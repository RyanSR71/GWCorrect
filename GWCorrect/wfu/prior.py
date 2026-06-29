import numpy as np
import bilby
import random
import time
import sys
import scipy
import lal
from .utils import A_ASD_solutions, gaussain_parameters_from_A_ASD_solutions
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
