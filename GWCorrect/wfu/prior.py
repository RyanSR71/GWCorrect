import numpy as np
import bilby
import random
import time
import sys
import scipy
import lal
from .utils import A_ASD_solutions, TFDG, EHG, xi_0_upper_bound, delta_xi_tilde_lower_bound
from bilby.core import utils
from bilby.core.series import CoupledTimeAndFrequencySeries
from bilby.core.utils import PropertyAccessor
from bilby.gw.conversion import convert_to_lal_binary_neutron_star_parameters



def xi_priors(waveform_generator,prior,psd_data,n,**kwargs):
    '''
    Generates xi_0 and delta_xi_tilde priors from a BBH/BNS/NSBH prior and adds them to the original prior.

    Parameters
    ==================
    waveform_generator: bilby.gw.WaveformGenerator
        bilby waveform generator object
    prior: bilby.core.prior.PriorDict
        bilby prior dictionary
    psd_data: numpy.ndarray
        array of power spectral density data; first column needs to be the frequency points and the second column needs to be the data
    n: int
        number of frequency nodes
    xi_0_latex_label: string, optional
        latex label for xi_0
        default: r'$xi_0$'
    delta_xi_tilde_latex_label: string, optional
        latex_label for delta_xi_tilde
        default: r'$delta tilde xi$'
    xi_min: float, optional
        lower bound on the dimensionless frequency band
        default: 0.018
    xi_max: float, optional
        upper bound on the dimensionless frequency band
        default: 1/pi
    samples: int, optional
        number of draws of amplitude to take to generate the priors
        default: 1000
    '''
    xi_0_latex_label = kwargs.get('xi_0_latex_label',r'$\xi_0$')
    delta_xi_tilde_latex_label = kwargs.get('delta_xi_tilde_latex_label',r'$\delta\tilde\xi$')
    xi_min = kwargs.get('xi_min',0.018)
    xi_max = kwargs.get('xi_max',1/np.pi)
    samples = kwargs.get('samples',1000)
    
    lower_xis, upper_xis = A_ASD_solutions(waveform_generator,psd_data,prior,samples,xi_min,xi_max,'Generating Priors')
    
    mu_1,sigma_1 = scipy.stats.norm.fit(lower_xis)
    mu_2,sigma_2 = scipy.stats.norm.fit(upper_xis)
    
    delta_xi_tildes = (np.array(upper_xis)-np.array(lower_xis))/(xi_max-np.array(lower_xis))
    
    mu_3,sigma_3 = scipy.stats.norm.fit(delta_xi_tildes)

    if mu_1 > xi_min:
        prior['xi_0'] = TFDG(name='xi_0',latex_label=xi_0_latex_label,
                            mu_1=mu_1,mu_2=mu_2,sigma_1=sigma_1,sigma_2=sigma_2,
                            minimum=xi_min,maximum=xi_max)
    else:
        prior['xi_0'] = EHG(name='xi_0',latex_label=xi_0_latex_label,
                            mu=mu_2,sigma=sigma_2,
                            minimum=xi_min,maximum=xi_max)
    
    prior['delta_xi_tilde'] = EHG(name='delta_xi_tilde',latex_label=delta_xi_tilde_latex_label,
                                  mu=mu_3,sigma=sigma_3,maximum=1,
                                  minimum=0)
    
    return prior



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



def conversion(parameters,xi_max,n):
    '''
    Conversion function to generate the total mass and frequency node constraint from a set of parameters;
    Necessary for the use of the constraint priors
    
    Parameters
    ==================
    parameters: dict
        dictionary of binary black hole parameters
    xi_max: float
        absolute upper bound on dimensionless frequency
    n: int
        number of frequency nodes (excluding the zeroth)
    
    Returns
    ==================
    parameters: dict
        input parameters, but with the constraint parameters added
    '''
    total_mass = bilby.gw.conversion.generate_mass_parameters(parameters)['total_mass']
    
    delta_t_01 = 1/((203025.4467280836/total_mass)*parameters['xi_0']*((1+((xi_max-parameters['xi_0'])/parameters['xi_0'])*parameters['delta_xi_tilde'])**(1/n)-1))

    parameters['total_mass'] = total_mass
    parameters['delta_t_01'] = delta_t_01
    
    return parameters
