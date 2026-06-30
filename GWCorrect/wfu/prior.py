import numpy as np
import bilby
import random
import time
import sys
import scipy
import lal
from .utils import A_ASD_solutions, gaussian_parameters_from_A_ASD_solutions, delta_t_01, epsilon_alpha
from bilby.core import utils
from bilby.core.series import CoupledTimeAndFrequencySeries
from bilby.core.utils import PropertyAccessor
from bilby.gw.conversion import convert_to_lal_binary_neutron_star_parameters



def TotalMassConstraint(*,name,minimum_frequency,maximum_frequency,**kwargs):
    '''
    Generates a bilby prior to constrain the total mass
    
    Parameters
    ===================
    name: string
        name of the prior
    minimum_frequency: float
        lower bound on the frequency band
    maximum_frequency: float
        upper bound of the frequency band
    unit: string, optional
        unit of the parameter
        default: r'$mathrm{M}_odot$' (solar mass)
    latex_label: string, optional
        label for the parameter in LaTeX
        default: r'$M$'
    boundary: string, optional
        boundary condition type for the prior
        default: None
    xi_min: float, optional
        lower bound on the waveform uncertainty correction in dimensionless frequency
        default: 0.018
    xi_max: float, optional
        upper bound on the waveform uncertainty correction in dimensionless frequency
        default: 1/pi, 0.318...

    Returns
    ==================
    total_mass_prior: bilby.core.prior.base.Constraint
        bilby constraint prior object for the total mass
    '''
    unit = kwargs.get('unit',r'$\mathrm{M}_\odot$')
    latex_label = kwargs.get('latex_label',r'$M$')
    boundary = kwargs.get('boundary',None)
    xi_min = kwargs.get('xi_min',0.018)
    xi_max = kwargs.get('xi_max',1/np.pi)
    
    total_mass_prior = bilby.core.prior.Constraint(name=name,latex_label=latex_label,minimum=xi_max*203025.4467280836/maximum_frequency, maximum=xi_min*203025.4467280836/minimum_frequency, unit=unit)
    
    return total_mass_prior



def conversion(input_parameters,**kwargs):
    '''
    Conversion function needed to use the total_mass constraint prior and the delta_t_01 constraint prior.

    Parameters
    ==================
    input_parameters: dict
        dictionary of binary black hole source parameters and waveform correction parameters
    n: int, optional
        number of frequency nodes excluding 0; if given, delta_t_01 will be calculated
        default: None
    xi_max: float, optional
        upper bound on the dimensionless frequency band
        default: 1/np.pi

    Returns
    ==================
    parameters: dict
        input parameters plus total mass and delta_t_01 (if n is not None)
    '''
    n = kwargs.get('n',None)
    xi_max = kwargs.get('xi_max',1/np.pi)
    
    parameters = input_parameters.copy()
    parameters['total_mass'] = bilby.gw.conversion.generate_mass_parameters(parameters)['total_mass']
    if n is not None:
        parameters['delta_t_01'] = delta_t_01(parameters['xi_0'],parameters['delta_xi_tilde'],parameters['total_mass'],n,xi_max)
    return parameters



def amplitude_conversion(input_parameters,n,xi_max,sigma_dA_spline):
    '''
    Conversion function for the amplitude error constraint parameter.

    Parameters
    ==================
    input_parameters: dict
        dictionary of binary black hole source parameters and waveform correction parameters
    n: int
        number of frequency nodes excluding 0; if given, delta_t_01 will be calculated
    xi_max: float
        upper bound on the dimensionless frequency band
    sigma_dA_spline: scipy.interpolate._cubic.CubicSpline
        scipy cubic spline object encoding the standard deviation of the dA prior

    Returns
    ==================
    parameters: dict
        input parameters plus epsilon_alpha
    '''
    parameters=input_parameters.copy()
    dAs = [parameters[f'dA_{k}'] for k in range(0,n+1)]
    parameters['epsilon_alpha'] = epsilon_alpha(dAs,parameters['xi_0'],parameters['delta_xi_tilde'],xi_max,sigma_dA_spline)
    return parameters



def frequency_node_prior_parameters(waveform_generator,asd_data,prior,**kwargs):
    '''
    Generates gaussian parameters for the xi_0 and delta_xi_tilde truncated gaussian priors by comparing waveforms to detector noise.

    Parameters
    ==================
    waveform_generator: bilby.gw.waveform_generator.WaveformGenerator
        bilby waveform generator object
    asd_data: numpy.ndarray
        array of amplitude spectral density data; assumes standard LIGO formatting
    prior: bilby.core.prior.PriorDict
        bilby prior object to draw BBH parameters from
    samples: int, optional
        number of waveforms to draw
        default: 10000
    xi_max: float, optional
        upper bound on the dimensionless frequency grid
        default: 1/pi, 0.318...

    Returns
    ==================
    mu_down: float
        mean of the xi_0 prior
    sigma_down: float
        standard deviation of the xi_0 prior
    mu: float
        mean of the delta_xi_tilde prior
    sigma: float
        standard deviation of the delta_xi_tilde prior
    '''
    samples = kwargs.get('samples',10000)
    xi_max = kwargs.get('xi_max',1/np.pi)
    
    lower_xis, upper_xis = A_ASD_solutions(waveform_generator,asd_data,prior,samples,'Generating Distributions')
    mu_down, sigma_down, _, _, mu, sigma = gaussian_parameters_from_A_ASD_solutions(lower_xis,upper_xis,xi_max)
    
    return mu_down, sigma_down, mu, sigma



class TFDG(bilby.core.prior.Prior):
    '''
    Truncated Flattened Double Gaussian Prior: a double Gaussian prior that is uniform between the two means; alternate prior for xi_0 if mu_down is greater than xi_min
        
    Parameters
    ==================
    mu_1: float
        mean for the first Gaussian
    sigma_1: float
        standard deviation for the first Gaussian
    mu_2: float
        mean for the second Gaussian
    sigma_2: float
        standard deviation for the second Gaussian
    minimum: float
        minimum value
    maximum: float
        maximum value

    Returns
    ==================
    TFDG prior object that is compatible with bilby.
    '''
    def __init__(self,mu_1,mu_2,sigma_1,sigma_2,minimum,maximum,name=None, latex_label=None):
        super(TFDG, self).__init__(
            name=name,latex_label=latex_label,minimum=minimum,maximum=maximum
        )
        self.mu_1 = float(mu_1)
        self.mu_2 = float(mu_2)
        self.sigma_1 = float(sigma_1)
        self.sigma_2 = float(sigma_2)
        
        
    def prob(self, val):
        in_region_1 = (val >= self.minimum) & (val <= self.mu_1)
        in_region_2 = (val > self.mu_1) & (val < self.mu_2)
        in_region_3 = (val >= self.mu_2) & (val <= self.maximum)
        N = (-np.sqrt(np.pi/2)*self.sigma_1*scipy.special.erf((self.minimum-self.mu_1)/(np.sqrt(2)*self.sigma_1))+np.sqrt(np.pi/2)*self.sigma_2*scipy.special.erf((self.maximum-self.mu_2)/(np.sqrt(2)*self.sigma_2))+self.mu_2-self.mu_1)**-1
        draw = N*np.exp(-0.5*((val-self.mu_1)/self.sigma_1)**2)*in_region_1+N*in_region_2+N*np.exp(-0.5*((val-self.mu_2)/self.sigma_2)**2)*in_region_3
        return draw
    
    
    def rescale(self, val):
        N = (-np.sqrt(np.pi/2)*self.sigma_1*scipy.special.erf((self.minimum-self.mu_1)/(np.sqrt(2)*self.sigma_1))+np.sqrt(np.pi/2)*self.sigma_2*scipy.special.erf((self.maximum-self.mu_2)/(np.sqrt(2)*self.sigma_2))+self.mu_2-self.mu_1)**-1
        A_1 = np.sqrt(np.pi/2)*N*self.sigma_1*scipy.special.erf((self.mu_1-self.minimum)/(np.sqrt(2)*self.sigma_1))
        A_2 = N*(self.mu_2-self.mu_1)
        
        if hasattr(val, "__len__"):
            draw = []
            for v in val:
        
                in_region_1 = (v >= 0) & (v <= A_1)
                in_region_2 = (v > A_1) & (v <= A_1+A_2)
                in_region_3 = (v >= A_1+A_2) & (v <= 1)

                if in_region_1:
                    draw.append(np.sqrt(2)*self.sigma_1*scipy.special.erfinv((np.sqrt(2*np.pi)*v-np.pi*np.sqrt(2/np.pi)*A_1)/(np.pi*N*self.sigma_1))+self.mu_1)
                elif in_region_2:
                    draw.append((N*self.mu_1+v-A_1)/N)
                elif in_region_3:
                    draw.append(self.mu_2-np.sqrt(2)*self.sigma_2*scipy.special.erfinv((np.sqrt(2*np.pi)*A_1+np.sqrt(2*np.pi)*A_2-np.sqrt(2*np.pi)*v)/(np.pi*N*self.sigma_2)))
                else:
                    raise Exception('Draw Failed!')
            return np.array(draw)
        
        else:
            in_region_1 = (val >= 0) & (val <= A_1)
            in_region_2 = (val > A_1) & (val <= A_1+A_2)
            in_region_3 = (val >= A_1+A_2) & (val <= 1)

            if in_region_1:
                draw = np.sqrt(2)*self.sigma_1*scipy.special.erfinv((np.sqrt(2*np.pi)*val-np.pi*np.sqrt(2/np.pi)*A_1)/(np.pi*N*self.sigma_1))+self.mu_1
            elif in_region_2:
                draw = (N*self.mu_1+val-A_1)/N
            elif in_region_3:
                draw = self.mu_2-np.sqrt(2)*self.sigma_2*scipy.special.erfinv((np.sqrt(2*np.pi)*A_1+np.sqrt(2*np.pi)*A_2-np.sqrt(2*np.pi)*val)/(np.pi*N*self.sigma_2))
            else:
                raise Exception('Draw Failed!')
            return draw



class LEHG(bilby.core.prior.Prior):
    '''
    Left Extended Half Gaussian Prior: a truncated Gaussian prior that is uniform between the minimum and the mean; alternate prior for xi_0 if mu_down is less than xi_min; alternate prior for delta_xi_tilde
        
    Parameters
    ==================
    mu: float
        mean
    sigma: float
        standard deviation
    minimum: float
        minimum value
    maximum: float
        maximum value

    Returns
    ==================
    LEHG prior object that is compatible with bilby.
    '''
    def __init__(self,mu,sigma,minimum,maximum,name=None, latex_label=None):
        super(LEHG, self).__init__(
            name=name,latex_label=latex_label,minimum=minimum,maximum=maximum
        )
        self.mu = float(mu)
        self.sigma = float(sigma)        
        
    def prob(self, val):
        in_region_1 = (val >= self.minimum) & (val <= self.mu)
        in_region_2 = (val > self.mu) & (val <= self.maximum)
        N = (np.sqrt(np.pi/2)*self.sigma*scipy.special.erf((self.maximum-self.mu)/(np.sqrt(2)*self.sigma))+self.mu-self.minimum)**-1
        draw = N*in_region_1+N*np.exp(-0.5*((val-self.mu)/self.sigma)**2)*in_region_2
        return draw
            
    
    def rescale(self, val):
        N = (np.sqrt(np.pi/2)*self.sigma*scipy.special.erf((self.maximum-self.mu)/(np.sqrt(2)*self.sigma))+self.mu-self.minimum)**-1
        A = N*(self.mu-self.minimum)
        
        if hasattr(val, "__len__"):
            draw = []
            for v in val:
                in_region_1 = (v >= 0) & (v < A)
                in_region_2 = (v >= A) & (v <= 1)

                if in_region_1:
                    draw.append((N*self.minimum+v)/N)
                elif in_region_2:
                    draw.append(self.mu-np.sqrt(2)*self.sigma*scipy.special.erfinv((np.sqrt(2*np.pi)*A-np.sqrt(2*np.pi)*v)/(np.pi*N*self.sigma)))
                else:
                    raise Exception('Draw Failed!')
            return np.array(draw)
        
        else:
            in_region_1 = (val >= 0) & (val < A)
            in_region_2 = (val >= A) & (val <= 1)

            if in_region_1:
                draw = ((N*self.minimum+val)/N)
            elif in_region_2:
                draw = (self.mu-np.sqrt(2)*self.sigma*scipy.special.erfinv((np.sqrt(2*np.pi)*A-np.sqrt(2*np.pi)*val)/(np.pi*N*self.sigma)))
            else:
                raise Exception('Draw Failed!')
            return draw



class REHG(bilby.core.prior.Prior):
    '''
    Right Extended Half Gaussian Prior: a truncated Gaussian prior that is uniform between the mean and the maximum; alternate prior for delta_xi_tilde
        
    Parameters
    ==================
    mu: float
        mean
    sigma: float
        standard deviation
    minimum: float
        minimum value
    maximum: float
        maximum value

    Returns
    ==================
    REHG prior object that is compatible with bilby.
    '''
    def __init__(self,mu,sigma,minimum,maximum,name=None, latex_label=None):
        super(REHG, self).__init__(
            name=name,latex_label=latex_label,minimum=minimum,maximum=maximum
        )
        self.mu = float(mu)
        self.sigma = float(sigma)        
        
    def prob(self, val):
        in_region_1 = (val >= self.minimum) & (val <= self.mu)
        in_region_2 = (val > self.mu) & (val <= self.maximum)
        N = (np.sqrt(np.pi/2)*self.sigma*scipy.special.erf((self.mu-self.minimum)/(np.sqrt(2)*self.sigma))+self.maximum-self.mu)**-1
        draw = N*in_region_2+N*np.exp(-0.5*((self.mu-val)/self.sigma)**2)*in_region_1
        return draw
            
    
    def rescale(self, val):
        N = (np.sqrt(np.pi/2)*self.sigma*scipy.special.erf((self.mu-self.minimum)/(np.sqrt(2)*self.sigma))+self.maximum-self.mu)**-1
        A = N*self.sigma*np.sqrt(np.pi/2)*scipy.special.erf((self.mu-self.minimum)/(np.sqrt(2)*self.sigma))
        
        if hasattr(val, "__len__"):
            draw = []
            for v in val:
                in_region_1 = (v >= 0) & (v < A)
                in_region_2 = (v >= A) & (v <= 1)

                if in_region_1:
                    draw.append(self.mu-np.sqrt(2)*self.sigma*scipy.special.erfinv((np.pi*N*self.sigma*scipy.special.erf((self.mu-self.minimum)/(np.sqrt(2)*self.sigma))-np.sqrt(2*np.pi)*v)/(np.pi*N*self.sigma)))
                elif in_region_2:
                    draw.append(self.mu+(v-A)/N)
                else:
                    raise Exception('Draw Failed!')
            return np.array(draw)
        
        else:
            in_region_1 = (val >= 0) & (val < A)
            in_region_2 = (val >= A) & (val <= 1)

            if in_region_1:
                draw = self.mu-np.sqrt(2)*self.sigma*scipy.special.erfinv((np.pi*N*self.sigma*scipy.special.erf((self.mu-self.minimum)/(np.sqrt(2)*self.sigma))-np.sqrt(2*np.pi)*val)/(np.pi*N*self.sigma))
            elif in_region_2:
                draw = self.mu+(val-A)/N
            else:
                raise Exception('Draw Failed!')
            return draw
