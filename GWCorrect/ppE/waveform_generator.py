import numpy as np
from numpy.linalg import inv
import bilby
import random
import time as tm
import sys
import scipy
import lal
import math
from .utils import inversion_function, beta_from_beta_tilde_wrapped, apply_ppe_correction

def ppECorrectionModel(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase,
        beta_tilde, delta_epsilon_tilde, b, **kwargs):
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
    kwargs: dict
        Optional keyword arguments
        Supported arguments:

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
    Mfreq_IM = kwargs.get('Mfreq_IM',0.018)
    ratio_f_MR_to_f_RD = kwargs.get('ratio_f_MR_to_f_RD',0.75)
    aligned = kwargs.get('aligned',True)
    
    waveform_arguments = dict(
        waveform_approximant=waveform_approximant, reference_frequency=reference_frequency,
        minimum_frequency=minimum_frequency, maximum_frequency=maximum_frequency,
        catch_waveform_errors=catch_waveform_errors, pn_spin_order=pn_spin_order, pn_tidal_order=pn_tidal_order,
        pn_phase_order=pn_phase_order, pn_amplitude_order=pn_amplitude_order,
    )
    
    unmodified_strain = bilby.gw.source._base_lal_cbc_fd_waveform(
        frequency_array=frequency_array, mass_1=mass_1, mass_2=mass_2,
        luminosity_distance=luminosity_distance, theta_jn=theta_jn, phase=phase,
        a_1=a_1, a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2, phi_12=phi_12,
        phi_jl=phi_jl, **waveform_arguments)
    
    beta = inversion_function(0,3,b)*beta_from_beta_tilde_wrapped(beta_tilde,minimum_frequency,1/np.pi,b,0.018,mass_1+mass_2)
    delta_epsilon = inversion_function(0,3,b)*beta_from_beta_tilde_wrapped(delta_epsilon_tilde,minimum_frequency,1/np.pi,b,0.018,mass_1+mass_2)
    
    model_strain = dict()
    model_strain['plus'] = apply_ppe_correction(unmodified_strain['plus'],frequency_array,mass_1+mass_2,beta,b,delta_epsilon,
                                                Mfreq_IM=Mfreq_IM,ratio_f_MR_to_f_RD=ratio_f_MR_to_f_RD,aligned=aligned,window_size=41)
    model_strain['cross'] = apply_ppe_correction(unmodified_strain['cross'],frequency_array,mass_1+mass_2,beta,b,delta_epsilon,
                                                 Mfreq_IM=Mfreq_IM,ratio_f_MR_to_f_RD=ratio_f_MR_to_f_RD,aligned=aligned,window_size=41)
    
    return model_strain



def MultippECorrectionModel(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase,
        beta_tildes, delta_epsilon_tildes, bs, **kwargs):
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
    kwargs: dict
        Optional keyword arguments
        Supported arguments:

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
    Mfreq_IM = kwargs.get('Mfreq_IM',0.018)
    ratio_f_MR_to_f_RD = kwargs.get('ratio_f_MR_to_f_RD',0.75)
    aligned = kwargs.get('aligned',True)
    
    waveform_arguments = dict(
        waveform_approximant=waveform_approximant, reference_frequency=reference_frequency,
        minimum_frequency=minimum_frequency, maximum_frequency=maximum_frequency,
        catch_waveform_errors=catch_waveform_errors, pn_spin_order=pn_spin_order, pn_tidal_order=pn_tidal_order,
        pn_phase_order=pn_phase_order, pn_amplitude_order=pn_amplitude_order,
    )
    
    unmodified_strain = bilby.gw.source._base_lal_cbc_fd_waveform(
        frequency_array=frequency_array, mass_1=mass_1, mass_2=mass_2,
        luminosity_distance=luminosity_distance, theta_jn=theta_jn, phase=phase,
        a_1=a_1, a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2, phi_12=phi_12,
        phi_jl=phi_jl, **waveform_arguments)

    model_strain = dict()
    model_strain['plus'] = unmodified_strain['plus']
    model_strain['cross'] = unmodified_strain['cross']

    n = len(bs)
    for i in range(n):

        b = bs[i]
        beta_tilde = beta_tildes[i]
        delta_epsilon_tilde = delta_epsilon_tildes[i]
        
        beta = inversion_function(0,3,b)*beta_from_beta_tilde_wrapped(beta_tilde,minimum_frequency,1/np.pi,b,0.018,mass_1+mass_2)
        delta_epsilon = inversion_function(0,3,b)*beta_from_beta_tilde_wrapped(delta_epsilon_tilde,minimum_frequency,1/np.pi,b,0.018,mass_1+mass_2)
    
        model_strain['plus'] = apply_ppe_correction(model_strain['plus'],frequency_array,mass_1+mass_2,beta,b,delta_epsilon,
                                                    Mfreq_IM=Mfreq_IM,ratio_f_MR_to_f_RD=ratio_f_MR_to_f_RD,aligned=aligned,window_size=41)
        model_strain['cross'] = apply_ppe_correction(model_strain['cross'],frequency_array,mass_1+mass_2,beta,b,delta_epsilon,
                                                     Mfreq_IM=Mfreq_IM,ratio_f_MR_to_f_RD=ratio_f_MR_to_f_RD,aligned=aligned,window_size=41)
    
    return model_strain



def binary_black_hole_ppE_conversion(parameters):
    """
    Convert parameters we have into parameters we need.

    This is defined by the parameters of bilby.source.lal_binary_black_hole()


    Mass: mass_1, mass_2
    Spin: a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl
    Extrinsic: luminosity_distance, theta_jn, phase, ra, dec, geocent_time, psi

    This involves popping a lot of things from parameters.
    The keys in added_keys should be popped after evaluating the waveform.

    Parameters
    ==========
    parameters: dict
        dictionary of parameter values to convert into the required parameters

    Returns
    =======
    converted_parameters: dict
        dict of the required parameters
    added_keys: list
        keys which are added to parameters during function call
    """

    converted_parameters = parameters.copy()
    original_keys = list(converted_parameters.keys())
    if 'luminosity_distance' not in original_keys:
        if 'redshift' in converted_parameters.keys():
            converted_parameters['luminosity_distance'] = \
                redshift_to_luminosity_distance(parameters['redshift'])
        elif 'comoving_distance' in converted_parameters.keys():
            converted_parameters['luminosity_distance'] = \
                comoving_distance_to_luminosity_distance(
                    parameters['comoving_distance'])

    for key in original_keys:
        if key[-7:] == '_source':
            if 'redshift' not in converted_parameters.keys():
                converted_parameters['redshift'] =\
                    luminosity_distance_to_redshift(
                        parameters['luminosity_distance'])
            converted_parameters[key[:-7]] = converted_parameters[key] * (
                1 + converted_parameters['redshift'])

    # we do not require the component masses be added if no mass parameters are present
    converted_parameters = bilby.gw.conversion.generate_component_masses(converted_parameters, require_add=False)

    for idx in ['1', '2']:
        key = 'chi_{}'.format(idx)
        if key in original_keys:
            if "chi_{}_in_plane".format(idx) in original_keys:
                converted_parameters["a_{}".format(idx)] = (
                    converted_parameters[f"chi_{idx}"] ** 2
                    + converted_parameters[f"chi_{idx}_in_plane"] ** 2
                ) ** 0.5
                converted_parameters[f"cos_tilt_{idx}"] = (
                    converted_parameters[f"chi_{idx}"]
                    / converted_parameters[f"a_{idx}"]
                )
            elif "a_{}".format(idx) not in original_keys:
                converted_parameters['a_{}'.format(idx)] = abs(
                    converted_parameters[key])
                converted_parameters['cos_tilt_{}'.format(idx)] = \
                    np.sign(converted_parameters[key])
            else:
                with np.errstate(invalid="raise"):
                    try:
                        converted_parameters[f"cos_tilt_{idx}"] = (
                            converted_parameters[key] / converted_parameters[f"a_{idx}"]
                        )
                    except (FloatingPointError, ZeroDivisionError):
                        logger.debug(
                            "Error in conversion to spherical spin tilt. "
                            "This is often due to the spin parameters being zero. "
                            f"Setting cos_tilt_{idx} = 1."
                        )
                        converted_parameters[f"cos_tilt_{idx}"] = 1.0

    for key in ["phi_jl", "phi_12"]:
        if key not in converted_parameters:
            converted_parameters[key] = 0.0

    for angle in ['tilt_1', 'tilt_2', 'theta_jn']:
        cos_angle = str('cos_' + angle)
        if cos_angle in converted_parameters.keys():
            with np.errstate(invalid="ignore"):
                converted_parameters[angle] = np.arccos(converted_parameters[cos_angle])

    if "delta_phase" in original_keys:
        with np.errstate(invalid="ignore"):
            converted_parameters["phase"] = np.mod(
                converted_parameters["delta_phase"]
                - np.sign(np.cos(converted_parameters["theta_jn"]))
                * converted_parameters["psi"],
                2 * np.pi)
    
    beta_tilde_keys = [key for key in parameters if 'beta_tilde_' in key]
    delta_epsilon_tilde_keys = [key for key in parameters if 'delta_epsilon_tilde_' in key]
    b_keys = [key for key in parameters if 'b_' in key]
    
    if len(beta_tilde_keys) != len(delta_epsilon_tilde_keys) != len(b_keys):
        raise Exception(f'Unequal numbers of ppE parameters given! beta_tildes: {len(beta_tilde_keys)} delta_epsilon_tildes: {len(delta_epsilon_tilde_keys)} bs: {len(b_keys)}')

    n = len(b_keys)
    
    beta_tildes = [parameters[beta_tilde_keys[i]] for i in range(n)]
    delta_epsilon_tildes = [parameters[delta_epsilon_tilde_keys[i]] for i in range(n)]
    bs = [parameters[b_keys[i]] for i in range(n)]

    converted_parameters['beta_tildes'] = beta_tildes
    converted_parameters['delta_epsilon_tildes'] = delta_epsilon_tildes
    converted_parameters['bs'] = bs
        
    added_keys = [key for key in converted_parameters.keys()
                  if key not in original_keys]

    return converted_parameters, added_keys
