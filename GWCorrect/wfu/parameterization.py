import numpy as np
import bilby
import scipy
import logging
import tqdm
from .utils import ProgressBar,convert_to_dimensionless_frequency



def fd_model_difference(hf1,hf2,**kwargs):
    '''
    Generates frequency domain waveform differences between two models hf1 and hf2
    
    Parameters
    ===================
    hf1: bilby.gw.waveform_generator.WaveformGenerator
        frequency domain waveform generator object
    hf2: bilby.gw.waveform_generator.WaveformGenerator
        frequency domain waveform generator object
    injection: dict, optional
        injection parameters if waveform generators do not have parameters in them
        default: None
    npoints: int, optional
        length of the desired frequency grid
        default: None; automatically set from frequency array
    polarization: string, optional
        polarization of the strain data (plus or cross)
        default: 'plus'
    psd_data: numpy.ndarray, optional
        array containing the psd data and their corresponding frequencies
        default: None
    correction_parameter: float, optional
        fraction of maximum amplitude to cut off the amplitude at
        default: 0.001
    ref_amplitude: numpy.ndarray, optional
        array of gravitational waveform amplitude
        default None
        
    Returns
    ==================
    frequency_grid: numpy.ndarray
        frequency grid that corresponds to the uncertainty arrays
    amplitude_difference: numpy.ndarray
        array of amplitude difference values
    phase_difference: numpy.ndarray
        array of phase difference values; if psd data is None, unaligned_phase_difference will be returned, aligned_phase_difference otherwise
    '''
    injection = kwargs.get('injection',None)
    npoints = kwargs.get('npoints',None)
    polarization = kwargs.get('polarization','plus')
    psd_data = kwargs.get('psd_data',None)
    correction_parameter = kwargs.get('correction_parameter',0.001)
    ref_amplitude = kwargs.get('ref_amplitude',None)

    minimum_frequency = hf1.waveform_arguments['minimum_frequency']
    maximum_frequency = hf1.waveform_arguments['maximum_frequency']
    
    bilby.core.utils.log.setup_logger(log_level=30)
    np.seterr(all='ignore')

    # adding injection parameters to waveform generators
    if injection is not None:
        hf1.frequency_domain_strain(parameters=injection)
        hf2.frequency_domain_strain(parameters=injection)
    
    # setting up frequency grid and frequency indexes
    start_index = np.argmin(np.abs(hf1.frequency_array - minimum_frequency))
    end_index = np.argmin(np.abs(hf1.frequency_array - maximum_frequency))
    if npoints is None:
        npoints = int(1+np.log(end_index/start_index)/np.log(1+1/start_index))
    elif npoints > int(1+np.log(end_index/start_index)/np.log(1+1/start_index)):
        raise ValueError(f'npoints value too large; maximum value is {int(1+np.log(end_index/start_index)/np.log(1+1/start_index))}')
    frequency_indexes = np.geomspace(start_index,end_index,npoints).astype(int)
    frequency_grid = hf1.frequency_array[frequency_indexes]

    # waveform amplitudes
    amplitude_1 = np.abs(hf1.frequency_domain_strain()[polarization][frequency_indexes])
    amplitude_2 = np.abs(hf2.frequency_domain_strain()[polarization][frequency_indexes])

    # waveform phases
    phase_1 = np.angle(hf1.frequency_domain_strain()[polarization][frequency_indexes])
    phase_2 = np.angle(hf2.frequency_domain_strain()[polarization][frequency_indexes])
                     
    amplitude_difference = (amplitude_2/amplitude_1) - 1
    unaligned_phase_difference = phase_2-phase_1

    # removing phase shifts of 2pi
    while any(value > 6 for value in [np.abs(unaligned_phase_difference[i+1]-unaligned_phase_difference[i]) for i in range(len(unaligned_phase_difference)-2)]):
        unaligned_phase_difference = np.unwrap(unaligned_phase_difference)

    # finding the position of the cuttoff/discontinuity correction frequency, f_disc
    final_index_1 = list(amplitude_1).index(min(amplitude_1, key=lambda x:np.abs(x-correction_parameter*np.max(amplitude_1))))
    final_index_2 = list(amplitude_2).index(min(amplitude_2, key=lambda x:np.abs(x-correction_parameter*np.max(amplitude_2))))
    final_index = min([final_index_1,final_index_2])
    
    # fitting a line to unaligned_phase_difference weighted by PSDs and subtracting off that line
    if psd_data is None:
        # loading in O4 asd data and converting to psd data
        asd_data = np.loadtxt('https://dcc.ligo.org/public/0165/T2000012/002/aligo_O4high.txt',comments='#')
        psd_data = asd_data.copy()
        psd_data[:,1] = psd_data[:,1]**2
    if ref_amplitude is None:
        ref_amplitude = np.abs(hf1.frequency_domain_strain()[polarization][frequency_indexes][0:final_index])
    ref_amplitude = np.interp(frequency_grid[0:final_index],maximum_frequency*np.linspace(0,1,len(ref_amplitude)),ref_amplitude)
    ref_sigma = np.interp(frequency_grid[0:final_index], psd_data[:,0],psd_data[:,1])
    align_weights = (ref_amplitude**2)/(ref_sigma)
    fit = np.polyfit(frequency_grid[0:final_index],unaligned_phase_difference[0:final_index],1,w=align_weights)
    unaligned_phase_difference_no_shifts = unaligned_phase_difference[0:final_index]-np.poly1d(fit)(frequency_grid[0:final_index])
    phase_difference = np.copy(unaligned_phase_difference)
    phase_difference[0:final_index] = unaligned_phase_difference_no_shifts
        
    # making the discontinuity correction to waveform model differences
    amplitude_difference[final_index:] = amplitude_difference[final_index-1]
    phase_difference[final_index:] = phase_difference[final_index-1]
    
    return frequency_grid,amplitude_difference,phase_difference



def parameterization(hf1,hf2,prior,nsamples,**kwargs):
    '''
    Generates samples of waveform differences between two models and parameterizes the data with Cubic splines.

    Parameters
    ==================
    hf1: bilby.gw.waveform_generator.WaveformGenerator
        frequency domain waveform generator object
    hf2: bilby.gw.waveform_generator.WaveformGenerator
        frequency domain waveform generator object
    prior: bilby.core.prior.dict.PriorDict
        bilby prior object
    nsamples: int
        number of draws of waveform uncertainty desired
    spline_resolution: float, optional
        fraction of frequency grid length to determine number of spline nodes; must be less than or equal to 1
        default: 0.5
    npoints: int, optional
        length of the desired frequency grid
        default: None; automatically set from frequency array
    psd_data: numpy.ndarray, optional
        array containing the psd data and their corresponding frequencies
        default: None
    correction_parameter: float, optional
        fraction of maximum amplitude to cut off the amplitude at
        default: 0.001
    ref_amplitude: numpy.ndarray, optional
        reference amplitude for residual phase calculation
        default: None
    polarization: string, optional
        polarization of the strain data (plus or cross)
        default: 'plus'

    Returns
    ==================
    parameterized_data: numpy.ndarray
        table containing the index, frequency_grid, dA_fit_parameters, dphi_fit_parameters, final_index, dA_final_point, dphi_final_point,
        and injection_parameters for each draw of waveform uncertainty

        frequency_grid: numpy.ndarray
            frequencies corresponding to the frequency parameters specified
        frequency_nodes: numpy.ndarray
            frequency nodes for the splines
        dA_parameters: numpy.ndarray
            amplitude difference spline parameters
        dphi_parameters: numpy.ndarray
            phase difference spline parameters
        injection: dictionary
            source parameters injected into the waveform generators
    '''
    npoints = kwargs.get('npoints',None)
    polarization = kwargs.get('polarization','plus')
    psd_data = kwargs.get('psd_data',None)
    correction_parameter = kwargs.get('correction_parameter',0.001)
    ref_amplitude = kwargs.get('ref_amplitude',None)
    spline_resolution = kwargs.get('spline_resolution',0.5)

    parameterization_data = np.zeros([nsamples,5],dtype=object)

    # setting the reference amplitude
    if ref_amplitude is None:
        injection = prior.sample()
        ref_amplitude = np.abs(hf1.frequency_domain_strain(parameters=injection)[polarization])

    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)
    log.addHandler(ProgressBar())
    
    for index in tqdm.tqdm(range(nsamples), desc = 'Generating Waveform Differences'):
        
        injection = prior.sample()
        
        # calculating waveform model differences
        frequency_grid,amplitude_difference,phase_difference = fd_model_difference(hf1,hf2,injection=injection,npoints=npoints,
                                                                                   polarization=polarization,psd_data=psd_data,
                                                                                   correction_parameter=correction_parameter,
                                                                                   ref_amplitude=ref_amplitude)

        spline_indexes = np.linspace(0,len(frequency_grid)-1,int(spline_resolution*len(frequency_grid))).astype(int)
        frequency_nodes = frequency_grid[spline_indexes]
        amplitude_parameters = amplitude_difference[spline_indexes]
        phase_parameters = phase_difference[spline_indexes]
        
        #constructing output matrix
        parameterization_data[index][0] = np.array(frequency_grid)
        parameterization_data[index][1] = np.array(frequency_nodes)
        parameterization_data[index][2] = np.array(amplitude_parameters)
        parameterization_data[index][3] = np.array(phase_parameters)
        parameterization_data[index][4] = injection
    
    return parameterization_data



def recovery_from_parameterization(parameterization_draw, **kwargs):
    '''
    Takes a row from a parameterization matrix and returns the frequency grid, amplitude difference, and phase difference

    Parameters
    ==================
    parameterization_draw: numpy.ndarray
        one row of a parameterization matrix
    dimensionless: bool, optional
        whether or not the output is returned in dimensionless frequency units
        Default: True
    xi_min: float, optional
        if dimensionless is True, this is the lower bound on the dimensionless frequency grid
        default: 0.001
    xi_max: float, optional
        if dimensionless is True, this is the upper bound on the dimensionless frequency grid
        default: 1/pi, 0.318...
    resolution: int, optional
        if dimensionless is True, this is the number of points in the dimensionless frequency grid
        default: 1000

    Returns:
    ==================
    frequency_grid: numpy.ndarray
        array of frequency points; if dimensionless=False, this will be the frequency grid stored in the parameterization file
    amplitude_difference: numpy.ndarray
        array of amplitude differences
    phase_difference: numpy.ndarray
        array of phase differences
    '''
    dimensionless = kwargs.get('dimensionless',True)
    xi_min = kwargs.get('xi_min',0.001)
    xi_max = kwargs.get('xi_max',1/np.pi)
    resolution = kwargs.get('resolution',1000)
    
    frequency_grid = parameterization_draw[0].copy()
    nodes = parameterization_draw[1].copy()
    dA_params = parameterization_draw[2].copy()
    dphi_params = parameterization_draw[3].copy()
    
    if dimensionless is True:
        total_mass = bilby.gw.conversion.generate_mass_parameters(parameterization_draw[4])['total_mass']
        total_frequency_grid = np.geomspace(xi_min,xi_max,resolution)
        frequency_grid = convert_to_dimensionless_frequency(frequency_grid,total_mass)
        nodes = convert_to_dimensionless_frequency(nodes,total_mass)
    else:
        total_frequency_grid = frequency_grid.copy()
        
    dA_spline = scipy.interpolate.CubicSpline(nodes,dA_params)
    amplitude_difference = np.interp(total_frequency_grid,frequency_grid,dA_spline(frequency_grid))
    
    dphi_spline = scipy.interpolate.CubicSpline(nodes,dphi_params)
    phase_difference = np.interp(total_frequency_grid,frequency_grid,dphi_spline(frequency_grid))
    
    return total_frequency_grid,amplitude_difference,phase_difference



def split_gaussian_parameters_from_parameterization(parameterization,**kwargs):
    '''
    Calculates the split-Gaussian parameters from a parameterization file for both the amplitude and phase differences.

    Parameters
    ==================
    parameterization: numpy.ndarray
        parameterization file from GWCorrect.wfu.parameterization.parameterization
    dimensionless: bool, optional
        whether or not the output is returned in dimensionless frequency units
        Default: True
    xi_min: float, optional
        if dimensionless is True, this is the lower bound on the dimensionless frequency grid
        default: 0.001
    xi_max: float, optional
        if dimensionless is True, this is the upper bound on the dimensionless frequency grid
        default: 1/pi, 0.318...
    resolution: int, optional
        if dimensionless is True, this is the number of points in the dimensionless frequency grid
        default: 1000

    Returns
    ==================
    dimensionless_frequency_grid: numpy.ndarray
        array of dimensionless frequency values corresponding to the output arrays
    dA_minus: numpy.ndarray
        standard deviation below the median amplitude difference
    median_dA: numpy.ndarray
        median amplitude difference values
    dA_plus: numpy.ndarray
        standard deviation above the median amplitude difference
    dphi_minus: numpy.ndarray
        standard deviation below the median phase difference
    median_dphi: numpy.ndarray
        median phase difference values
    dphi_plus: numpy.ndarray
        standard deviation above the median phase difference    
    '''
    dimensionless = kwargs.get('dimensionless',True)
    xi_min = kwargs.get('xi_min',0.001)
    xi_max = kwargs.get('xi_max',1/np.pi)
    resolution = kwargs.get('resolution',1000)
    
    amplitude_differences = []
    phase_differences = []
    for i in range(len(parameterization)):
        dimensionless_frequency_grid,amplitude_difference,phase_difference = recovery_from_parameterization(parameterization[i],dimensionless=dimensionless,
                                                                                                            xi_min=xi_min,xi_max=xi_max,resolution=resolution)
        
        amplitude_differences.append(amplitude_difference)
        phase_differences.append(phase_difference)
    
    dA_plus = []
    dA_minus = []
    median_dA = []
    
    dphi_plus = []
    dphi_minus = []
    median_dphi = []
    
    for i in range(len(amplitude_differences[0])):
        dA_data = []
        dphi_data = []
        for j in range(len(amplitude_differences)):
            dA_data.append(amplitude_differences[j][i])
            dphi_data.append(phase_differences[j][i])
        dA_data = np.array(dA_data)
        dphi_data = np.array(dphi_data)

        dA_15,dA_50,dA_84 = np.percentile(dA_data,np.array([15.865,50,84.135]))
        dphi_15,dphi_50,dphi_84 = np.percentile(dphi_data,np.array([15.865,50,84.135]))

        dA_plus.append(dA_84-dA_50)
        dA_minus.append(dA_50-dA_15)
        median_dA.append(dA_50)
        
        dphi_plus.append(dphi_84-dphi_50)
        dphi_minus.append(dphi_50-dphi_15)
        median_dphi.append(dphi_50)
        
    return dimensionless_frequency_grid, np.array(dA_minus), np.array(median_dA), np.array(dA_plus), np.array(dphi_minus), np.array(median_dphi), np.array(dphi_plus)



def sigma_splines_from_parameterization(parameterization,**kwargs):
    '''
    Generates splines encoding the total standard deviations for the amplitude and phase correction priors

    Parameters
    ==================
    parameterization: numpy.ndarray
        parameterization file from GWCorrect.wfu.parameterization.parameterization
    dimensionless: bool, optional
        whether or not the output is returned in dimensionless frequency units
        Default: True
    xi_min: float, optional
        if dimensionless is True, this is the lower bound on the dimensionless frequency grid
        default: 0.001
    xi_max: float, optional
        if dimensionless is True, this is the upper bound on the dimensionless frequency grid
        default: 1/pi, 0.318...
    resolution: int, optional
        if dimensionless is True, this is the number of points in the dimensionless frequency grid
        default: 1000

    Returns
    ==================
    sigma_A_spline: scipy.interpolate._cubic.CubicSpline
        scipy cubic spline object encoding the total standard deviation of the amplitude correction prior
    sigma_phi_spline: scipy.interpolate._cubic.CubicSpline
        scipy cubic spline object encoding the total standard deviation of the phase correction prior
    '''
    dimensionless = kwargs.get('dimensionless',True)
    xi_min = kwargs.get('xi_min',0.001)
    xi_max = kwargs.get('xi_max',1/np.pi)
    resolution = kwargs.get('resolution',1000)
    
    dimensionless_frequency_grid, dA_minus, median_dA, dA_plus, dphi_minus, median_dphi, dphi_plus = split_gaussian_parameters_from_parameterization(parameterization,
                                                                                                                                                     dimensionless=dimensionless,
                                                                                                                                                     xi_min=xi_min,xi_max=xi_max,
                                                                                                                                                     resolution=resolution)
        
    total_dA_variance = (1-2/np.pi)*(dA_plus-dA_minus)**2+dA_plus*dA_minus
    total_dphi_variance = (1-2/np.pi)*(dphi_plus-dphi_minus)**2+dphi_plus*dphi_minus
    
    sigma_A = np.sqrt(median_dA**2+total_dA_variance)
    sigma_phi = np.sqrt(median_dphi**2+total_dphi_variance)
    
    spline_nodes = dimensionless_frequency_grid.copy()
    
    sigma_A_spline = scipy.interpolate.CubicSpline(spline_nodes,sigma_A)
    sigma_phi_spline = scipy.interpolate.CubicSpline(spline_nodes,sigma_phi)
    
    return sigma_A_spline, sigma_phi_spline
