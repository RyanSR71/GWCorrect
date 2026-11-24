import numpy as np
import bilby
import logging
import tqdm
import pycbc
import pycbc.filter
from ..wfu.utils import ProgressBar
import matplotlib.pyplot as plt
from pesummary.gw.file.strain import StrainData
from pesummary.io import read

def match(waveform_1,waveform_2,injection):
    '''
    Computes the match between two waveforms

    Parameters
    ==================
    waveform_1: bilby.gw.waveform_generator.WaveformGenerator
        bilby waveform generator object
    waveform_2: bilby.gw.waveform_generator.WaveformGenerator
        bilby waveform generator object
    injection: dict
        dictionary of injection parameters to generate waveforms with the waveform generators

    Returns
    ==================
    match: float
        match between the two waveforms at the given injection
    '''
    delta_f = waveform_1.frequency_array[1]-waveform_1.frequency_array[0]
    strain_1 = pycbc.types.frequencyseries.FrequencySeries(waveform_1.frequency_domain_strain(parameters=injection)['plus'],delta_f=delta_f)
    strain_2 = pycbc.types.frequencyseries.FrequencySeries(waveform_2.frequency_domain_strain(parameters=injection)['plus'],delta_f=delta_f)
    match = pycbc.filter.matchedfilter.match(strain_1,strain_2)[0]
    return match



def match_plot(ppE_waveform_generator,GR_waveform_generator,injection,beta_tildes,delta_epsilon_tildes):
    '''
    Generates contour plot parameters for displaying the match between a ppE waveform and a GR waveform at different values of beta_tilde and delta_epsilon_tilde.

    Parameters
    ==================
    ppE_waveform_generator: bilby.gw.waveform_generator.WaveformGenerator
        bilby waveform generator object with GWCorrect.ppE.waveform_generator.ppECorrectionModel as the frequency domain source model
    GR_waveform_generator: bilby.gw.waveform_generator.WaveformGenerator
        default bilby waveform generator
    injection: dict
        dictionary of BBH injection parameters; include a value for the ppE PN order parameter, b
    beta_tildes: numpy.ndarray
        array of beta_tilde values to evaluate the matches over
    delta_epsilon_tildes: numpy.ndarray
        array of delta_epsilon_tilde values to evaluate the matches over

    Returns
    ==================
    X: numpy.ndarray
        matrix of beta_tildes to use with matplotlib.pyplot.contour
    Y: numpy.ndarray
        matrix of delta_epsilon_tildes to use with matplotlib.pyplot.contour
    match_grid: numpy.ndarray
        matrix of match values corresponding to the injected beta_tildes and delta_epsilon_tildes; used with matplotlib.pyplot.contour
    '''
    if len(beta_tildes) != len(delta_epsilon_tildes):
        raise Exception(f'beta_tilde array and delta_epsilon_tilde array do not have the same length! {len(beta_tildes)}, {len(delta_epsilon_tildes)}')
        
    resolution = len(beta_tildes)
    match_grid = np.zeros([resolution,resolution])
    
    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)
    log.addHandler(ProgressBar())
    with tqdm.tqdm(total=resolution**2, desc="Generating Contour Plot",position=0,leave=True) as pbar:
        for i in range(resolution):
            for j in range(resolution):

                injection['beta_tilde'] = beta_tildes[i]
                injection['delta_epsilon_tilde'] = delta_epsilon_tildes[j]

                match_value = match(ppE_waveform_generator,GR_waveform_generator,injection)

                match_grid[j,i] = match_value

                pbar.update(1)
                
    X,Y = np.meshgrid(beta_tildes,delta_epsilon_tildes)

    return X, Y, match_grid



def total_mass_conversion(parameters):
    '''
    Conversion function to generate the total mass from a set of parameters; to be used alongside the total mass prior
    
    Parameters
    ==================
    parameters: dict
        dictionary of binary black hole parameters
    
    Returns
    ==================
    parameters: dict
        input parameters, but with the total mass added
    '''
    parameters['total_mass'] = bilby.gw.conversion.generate_mass_parameters(parameters)['total_mass']
    
    return parameters



def TotalMassConstraintPPE(*,name,minimum_frequency,**kwargs):
    '''
    Generates a bilby prior to constrain the total mass
    
    Parameters
    ===================
    name: string
        name of the prior
    minimum_frequency: float
        lower bound on the frequency band
    unit: string, optional
        unit of the parameter
        default: r'$mathrm{M}_odot$' (solar mass)
    latex_label: string, optional
        label for the parameter in LaTeX
        default: r'$M$'
    boundary: string, optional
        boundary condition type for the prior
        default: None
    Mf_IM: float, optional
        end of the early inspiral regime in geometrized units
        default: 0.018

    Returns
    ==================
    total_mass_prior: bilby.core.prior.base.Constraint
        bilby constraint prior object for the total mass
    '''
    unit = kwargs.get('unit',r'$\mathrm{M}_\odot$')
    latex_label = kwargs.get('latex_label',r'$M$')
    boundary = kwargs.get('boundary',None)
    Mf_IM = kwargs.get('xi_low',0.018)
    
    total_mass_prior = bilby.core.prior.Constraint(name=name,latex_label=latex_label,minimum=0, maximum=Mf_IM*203025.4467280836/minimum_frequency, unit=unit)
    
    return total_mass_prior
