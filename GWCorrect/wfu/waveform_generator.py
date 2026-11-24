import numpy as np
import bilby
import scipy
import logging
from bilby.core.utils import logger
from .utils import amplitude_difference_from_parameters, phase_difference_from_parameters
from .conversion import dimensionless_frequency_nodes_from_delta_xis



def WaveformCorrectionModelBBH(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase,
        xi_0, delta_xi_tilde, dAs, dphis, **kwargs):
    """ A Binary Black Hole waveform model using lalsimulation

    Parameters
    ==========
    frequency_array: array_like
        The frequencies at which we want to calculate the strain
    mass_1: float
        The mass of the heavier object in solar masses
    mass_2: float
        The mass of the lighter object in solar masses
    luminosity_distance: float
        The luminosity distance in megaparsec
    a_1: float
        Dimensionless primary spin magnitude
    tilt_1: float
        Primary tilt angle
    phi_12: float
        Azimuthal angle between the two component spins
    a_2: float
        Dimensionless secondary spin magnitude
    tilt_2: float
        Secondary tilt angle
    phi_jl: float
        Azimuthal angle between the total binary angular momentum and the
        orbital angular momentum
    theta_jn: float
        Angle between the total binary angular momentum and the line of sight
    phase: float
        The phase at reference frequency or peak amplitude (depends on waveform)
    dimensionless_frequency_nodes: list
        list of dimensionless_frequency_nodes
    alpha_tildes: list
        list of rescaled amplitude correction parameters
    phi_tildes: list
        list of rescaled phase correction parameters
    
    kwargs: dict
        Optional keyword arguments
        Supported arguments:

        - sigma_A_spline:
          total standard deviation of the amplitude correction prior
        - sigma_phi_spline:
          total standard deviation of the phase correction prior
        - gamma
          smoothing parameter for the smooth interpolation functions; default: 0.025

        - waveform_approximant
        - reference_frequency
        - minimum_frequency
        - maximum_frequency
        - catch_waveform_errors
        - pn_spin_order
        - pn_tidal_order
        - pn_phase_order
        - pn_amplitude_order
        - mode_array:
          Activate a specific mode array and evaluate the model using those
          modes only.  e.g. waveform_arguments =
          dict(waveform_approximant='IMRPhenomHM', mode_array=[[2,2],[2,-2])
          returns the 22 and 2-2 modes only of IMRPhenomHM.  You can only
          specify modes that are included in that particular model.  e.g.
          waveform_arguments = dict(waveform_approximant='IMRPhenomHM',
          mode_array=[[2,2],[2,-2],[5,5],[5,-5]]) is not allowed because the
          55 modes are not included in this model.  Be aware that some models
          only take positive modes and return the positive and the negative
          mode together, while others need to call both.  e.g.
          waveform_arguments = dict(waveform_approximant='IMRPhenomHM',
          mode_array=[[2,2],[4,-4]]) returns the 22 and 2-2 of IMRPhenomHM.
          However, waveform_arguments =
          dict(waveform_approximant='IMRPhenomXHM', mode_array=[[2,2],[4,-4]])
          returns the 22 and 4-4 of IMRPhenomXHM.
        - lal_waveform_dictionary:
          A dictionary (lal.Dict) of arguments passed to the lalsimulation
          waveform generator. The arguments are specific to the waveform used.

    Returns
    =======
    dict: A dictionary with the plus and cross polarisation strain modes
    """
    waveform_approximant = kwargs.get('waveform_approximant','IMRPhenomPv2')
    reference_frequency = kwargs.get('reference_frequency',50.0)
    minimum_frequency = kwargs.get('minimum_frequency',20.0)
    maximum_frequency = kwargs.get('maximum_frequency',frequency_array[-1])
    catch_waveform_errors = kwargs.get('catch_waveform_errors',False)
    pn_spin_order = kwargs.get('pn_spin_order',-1)
    pn_tidal_order = kwargs.get('pn_tidal_order',-1)
    pn_phase_order = kwargs.get('pn_phase_order',-1)
    pn_amplitude_order = kwargs.get('pn_amplitude_order',0)
    sigma_A_spline = kwargs.get('sigma_A_spline',None)
    sigma_phi_spline = kwargs.get('sigma_phi_spline',None)
    gamma = kwargs.get('gamma',0.025)
    
    waveform_arguments = dict(
        waveform_approximant=waveform_approximant, reference_frequency=reference_frequency,
        minimum_frequency=minimum_frequency, maximum_frequency=maximum_frequency,
        catch_waveform_errors=catch_waveform_errors, pn_spin_order=pn_spin_order, pn_tidal_order=pn_tidal_order,
        pn_phase_order=pn_phase_order, pn_amplitude_order=pn_amplitude_order,
    )
    
    model_strain = bilby.gw.source._base_lal_cbc_fd_waveform(
        frequency_array=frequency_array, mass_1=mass_1, mass_2=mass_2,
        luminosity_distance=luminosity_distance, theta_jn=theta_jn, phase=phase,
        a_1=a_1, a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2, phi_12=phi_12,
        phi_jl=phi_jl, **waveform_arguments)

    if sigma_A_spline is not None:
        amplitude_difference = amplitude_difference_from_parameters(frequency_array,dimensionless_frequency_nodes,alpha_tildes,sigma_A_spline,mass_1,mass_2,gamma)
    else:
        amplitude_difference = 0

    if sigma_phi_spline is not None:
        phase_difference = phase_difference_from_parameters(frequency_array,dimensionless_frequency_nodes,phi_tildes,sigma_phi_spline,mass_1,mass_2,gamma)
    else:
        phase_difference = 0

    waveform_correction = (1+amplitude_difference)*np.exp(phase_difference*1j)
    
    model_strain['plus'] *= waveform_correction
    model_strain['cross'] *= waveform_correction
    
    return model_strain



def binary_black_hole_correction_conversion(parameters):
    '''
    Converts source and waveform correction parameters to use with WaveformCorrectionModelBBH

    Parameters
    ==================
    parameters: dict
        dictionary of input parameters

    Returns
    ==================
    converted_parameters: dict
        dictionary of converted parameters
    added_keys: list
        list of parameters added to the parameter dictionary
    '''
    converted_parameters,_ = bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters(parameters)

    delta_xi_keys = [key for key in parameters if 'delta_xi_' in key]
    alpha_tilde_keys = [key for key in parameters if 'alpha_tilde_' in key]
    phi_tilde_keys = [key for key in parameters if 'phi_tilde_' in key]

    if 'xi_0' and 'xi_n' in parameters.keys():
        n = len(delta_xi_keys) + 2 #add an extra 1 to account for xi_0
    else:
        n = 0

    if n > 2:
        delta_xis = [parameters[key] for key in delta_xi_keys]
        dimensionless_frequency_nodes = dimensionless_frequency_nodes_from_delta_xis(parameters['xi_0'],parameters['xi_n'],delta_xis)
    elif n == 2:
        dimensionless_frequency_nodes = [parameters['xi_0'],parameters['xi_n']]
    else:
        dimensionless_frequency_nodes = [0.01,1] #arbitrary choice

    converted_parameters['dimensionless_frequency_nodes'] = dimensionless_frequency_nodes

    if 'alpha_tilde_0' in alpha_tilde_keys:
        correction = 0
    elif 'phi_tilde_0' in phi_tilde_keys:
        correction = 0
    else:
        correction = 1
    
    if len(alpha_tilde_keys) == 0:
        correct_amplitude = False
    elif len(alpha_tilde_keys)+correction == n:
        alpha_tildes = [parameters[key] for key in alpha_tilde_keys]
        if not any(alpha_tildes):
            correct_amplitude = False
        else:
            correct_amplitude = True
    else:
        raise Exception(f'Amplitude correction does not have the same number of parameters as frequency nodes! alpha_tildes: {len(alpha_tilde_keys)} n: {n}')
    
    if len(phi_tilde_keys) == 0:
        correct_phase = False
    elif len(phi_tilde_keys)+correction == n:
        phi_tildes = [parameters[key] for key in phi_tilde_keys]
        if not any(phi_tildes):
            correct_phase = False
        else:
            correct_phase = True
    else:
        raise Exception(f'Phase correction does not have the same number of parameters as frequency nodes! phi_tildes: {len(phi_tilde_keys)} n: {n}')
    
    if correct_amplitude and correct_phase:
            converted_parameters['alpha_tildes'] = alpha_tildes
            converted_parameters['phi_tildes'] = phi_tildes
    
    if correct_amplitude and not correct_phase:
        converted_parameters['alpha_tildes'] = alpha_tildes
        converted_parameters['phi_tildes'] = list(np.zeros(len(alpha_tildes)))
        
    if correct_phase and not correct_amplitude:
        converted_parameters['phi_tildes'] = phi_tildes
        converted_parameters['alpha_tildes'] = list(np.zeros(len(phi_tildes)))
        
    if not (correct_phase or correct_amplitude):
        converted_parameters['alpha_tildes'] = [0,0]
        converted_parameters['phi_tildes'] = [0,0]
        
    added_keys = [key for key in converted_parameters.keys()]

    return converted_parameters, added_keys
