import numpy as np
import bilby
import scipy
import tqdm
import logging
from .conversion import dimensionless_frequency_nodes_from_xi_0_delta_xi_tilde,convert_to_dimensionless_frequency



class ProgressBar(logging.Handler):
    '''
    Progress bar utility function. To use, set this as your logger handler.
    '''
    def __init__(self, level=logging.NOTSET):
        super().__init__(level)

    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.tqdm.write(msg)
            self.flush()
        except Exception:
            self.handleError(record)



def rolling_average(data,gamma):
    '''
    Rolling average function.

    Parameters
    ==================
    data: numpy.ndarray
        input array to be averaged
    gamma: float
        smoothing parameter
    '''
    new_data = data.copy()
    r = int(gamma*len(data))
    if r != 0:
        lower_index = int(gamma*len(data))
        upper_index = int((1-gamma)*len(data))
        for i in range(lower_index,upper_index):
            new_data[i] = (1/(2*r))*np.sum(data[i-r:i+r])

    return new_data



def smooth_interpolation(full_grid,nodes,parameters,gamma):
    '''
    Smooth interpolation function.

    Parameters
    ==================
    full_grid: numpy.ndarray
        array of points over which the smooth interpolation function is generated
    nodes: numpy.ndarray
        array of nodes (x_i)
    parameters: numpy.ndarray
        array of parameters corresponding with the nodes (y_i)
    gamma: float
        smoothing parameter

    Returns
    ==================
    smooth_interpolation: numpy.ndarray
        smooth_interpolation array that corresponds to the full_grid input
    '''
    spline = scipy.interpolate.interp1d(nodes,parameters)(nodes)
    temp_grid = np.geomspace(full_grid[1],full_grid[-1],200)
    data = np.interp(temp_grid,nodes,spline)

    new_data = rolling_average(data,gamma)
    
    output = np.interp(full_grid,temp_grid,new_data)
    if np.abs(output[0]) > 0:
        output -= output[0]

    return output



def dh_SI(dimensionless_frequency_grid,dimensionless_frequency_nodes,dhs,sigma_dh_spline,**kwargs):
    '''
    Smooth interpolation function for the waveform correction parameters.

    Parameters
    ==================
    dimensionless_frequency_grid: numpy.ndarray
        array of mass-weighted dimensionless frequencies
    dimensionless_frequency_nodes: list or numpy.ndarray
        list of dimensionless frequency nodes
    dhs: list or numpy.ndarray
        rescaled waveform correction parmaters, either dA or dphi
    sigma_dh_spline: scipy.interpolate._cubic.CubicSpline
        scipy cubic spline object encoding the standard deviation of the dh prior, either dA or dphi
    gamma: float, optional
        smoothing parameters
        default: 0.025

    Returns
    ==================
    dh_smooth_interpolation: numpy.ndarray
        dA or dphi smooth interpolation array corresponding to the input dimensionless frequency grid
    '''
    gamma = kwargs.get('gamma',0.025)
    n = len(dhs)-1
    true_dhs = [dhs[k]*sigma_dh_spline(dimensionless_frequency_nodes[k]) for k in range(0,n+1)]
    dh_SI = smooth_interpolation(dimensionless_frequency_grid,dimensionless_frequency_nodes,true_dhs,gamma)
    return dh_SI



def dh_1D(dimensionless_frequency_grid,dimensionless_frequency_nodes,dhs,sigma_dh_spline):
    '''
    1D interpolation function for the waveform correction parameters.

    Parameters
    ==================
    dimensionless_frequency_grid: numpy.ndarray
        array of mass-weighted dimensionless frequencies
    dimensionless_frequency_nodes: list or numpy.ndarray
        list of dimensionless frequency nodes
    dhs: list or numpy.ndarray
        rescaled waveform correction parmaters, either dA or dphi
    sigma_dh_spline: scipy.interpolate._cubic.CubicSpline
        scipy cubic spline object encoding the standard deviation of the dh prior, either dA or dphi

    Returns
    ==================
    dh_1D_interpolation: numpy.ndarray
        dA or dphi 1D interpolation array corresponding to the input dimensionless frequency grid
    '''
    n = len(dhs)-1
    true_dhs = [dhs[k]*sigma_dh_spline(dimensionless_frequency_nodes[k]) for k in range(0,n+1)]
    dh_1D_object = scipy.interpolate.interp1d(dimensionless_frequency_nodes,true_dhs)
    return dh_1D_object(dimensionless_frequency_grid)



def remove_dphi_shifts(dimensionless_frequency_grid,injection,n,sigma_dphi_spline,asd_data,waveform_generator,**kwargs):
    '''
    Takes an injection of binary black hole and phase correction parameters and removes any time domain shifts from the phase deviation.

    Parameters
    ==================
    dimensionless_frequency_grid: numpy.ndarray
        array of dimensionless frequencies to evaluate the phase difference over
    injection: dict
        dictionary of binary black hole source parameters and phase correction parameters
    n: int
        number of frequency nodes excluding 0
    sigma_dphi_spline: scipy.interpolate._cubic.CubicSpline
        scipy cubic spline object encoding the standard deviation of the dphi prior
    asd_data: numpy.ndarray
        array of amplitude spectral density data; assumes this is in the standard LIGO format
    waveform_generator: bilby.gw.waveform_generator.WaveformGenerator
        bilby waveform generator object
    xi_max: float, optional
        upper bound on the dimensionless frequency band
        default: 1/pi, 0.318...
    gamma: float, optional
        smoothing parameter
        default: 0.025

    Returns
    ==================
    phase_difference_no_shifts: numpy.ndarray
        input phase difference from injection with overall time and phase shifts removed; same shape as input dimensionless frequency grid
    '''
    xi_max = kwargs.get('xi_max',1/np.pi)
    gamma = kwargs.get('gamma',0.025)
    
    dimensionless_frequency_nodes = dimensionless_frequency_nodes_from_xi_0_delta_xi_tilde(injection['xi_0'],injection['delta_xi_tilde'],3,xi_max)
    temp_dimensionless_frequency_grid = np.geomspace(dimensionless_frequency_nodes[0],dimensionless_frequency_nodes[-1],1000)
    
    dphis = [injection[f'dphi_{k}'] for k in range(0,n+1)]
    phase_difference_1D = dh_1D(temp_dimensionless_frequency_grid,dimensionless_frequency_nodes,dphis,sigma_dphi_spline)
    
    reference_amplitude = np.sqrt(waveform_generator.frequency_domain_strain(parameters=injection)['plus']**2+waveform_generator.frequency_domain_strain(parameters=injection)['cross']**2)
    total_mass = bilby.gw.conversion.generate_mass_parameters(injection)['total_mass']
    reference_amplitude_interp = np.interp(temp_dimensionless_frequency_grid,convert_to_dimensionless_frequency(waveform_generator.frequency_array,total_mass),reference_amplitude)
    psd_data_interp = np.interp(temp_dimensionless_frequency_grid,convert_to_dimensionless_frequency(asd_data[:,0],total_mass),asd_data[:,1]**2)
    weights = np.abs(reference_amplitude_interp**2 / psd_data_interp)
    fit = np.polyfit(temp_dimensionless_frequency_grid,phase_difference_1D,1,w=weights)
    
    phase_difference_no_shifts_1D = phase_difference_1D - np.poly1d(fit)(temp_dimensionless_frequency_grid)
    phase_difference_no_shifts_interp = np.interp(dimensionless_frequency_grid,temp_dimensionless_frequency_grid,phase_difference_no_shifts_1D)
    phase_difference_no_shifts = rolling_average(phase_difference_no_shifts_interp,gamma)
    
    return phase_difference_no_shifts



def BBH_waveform_correction(frequency_array,xi_0,delta_xi_tilde,dAs,dphis,sigma_dA_spline,sigma_dphi_spline,mass_1,mass_2,xi_max,gamma):
    '''
    Computes (1+dA)exp(idphi) from the frequency node and correction parameters.

    Parameters
    ==================
    frequency_array: numpy.ndarray
        frequency grid over which to generate the waveform correction
    xi_0: float
        starting frequency node parameter
    delta_xi_tilde: float
        frequency node spacing parameter
    dAs: list or numpy.ndarray
        list of dA (alpha_tilde) parameters
    dphis: list or numpy.ndarray
        list of dphi (varphi_tilde) parameters
    sigma_dA_spline: scipy.interpolate._cubic.CubicSpline
        scipy cubic spline object encoding the standard deviation of the dA prior
    sigma_dphi_spline: scipy.interpolate._cubic.CubicSpline
        scipy cubic spline object encoding the standard deviation of the dphi prior
    mass_1: float
        mass of the primary body in the binary
    mass_2: float
        mass of the secondary body in the binary
    xi_max: float
        upper bound on the dimensionless frequency grid
    gamma: float
        smoothing parameter for the smooth interpolation function

    Returns
    ==================
    waveform_correction: numpy.ndarray
        (1+dA)exp(idphi) array to multiply to frequency domain gravitational-wave strain; same shape as frequency_array input
    '''
    n = len(dAs)-1
    total_mass = mass_1+mass_2
    dimensionless_frequency_nodes = np.array([xi_0*(1+((xi_max-xi_0)/(xi_0))*delta_xi_tilde)**(k/n) for k in range(n+1)])
    dimensionless_frequency_array = np.linspace(dimensionless_frequency_nodes[0],dimensionless_frequency_nodes[-1],len(frequency_array))
    
    if sigma_dA_spline is not None:
        sigma_dA = sigma_dA_spline(dimensionless_frequency_nodes)
    else: 
        sigma_dA = np.ones(n+1)
    if sigma_dphi_spline is not None:
        sigma_dphi = sigma_dphi_spline(dimensionless_frequency_nodes)
    else:
        sigma_dphi = np.ones(n+1)
    
    amplitude_correction = smooth_interpolation(dimensionless_frequency_array,dimensionless_frequency_nodes,dAs*sigma_dA,gamma)
    amplitude_correction = np.interp(frequency_array,dimensionless_frequency_array*203025.4467280836/total_mass,amplitude_correction)
    
    phase_correction = smooth_interpolation(dimensionless_frequency_array,dimensionless_frequency_nodes,dphis*sigma_dphi,gamma)
    phase_correction = np.interp(frequency_array,dimensionless_frequency_array*203025.4467280836/total_mass,phase_correction)
    
    waveform_correction = (1+amplitude_correction)*np.exp(phase_correction*1j)
    
    return waveform_correction



def A_ASD_solutions(waveform_generator,asd_data,prior,samples,desc):
    '''
    Computes the distributions in dimensionless frequency where waveforms cross detector noise.

    Parameters
    ==================
    waveform_generator: bilby.gw.waveform_generator.WaveformGenerator
        bilby waveform generator object
    asd_data: numpy.ndarray
        array of amplitude spectral density data; assumes standard LIGO formatting
    prior: bilby.core.prior.PriorDict
        bilby prior object to draw BBH parameters from
    samples: int
        number of waveforms to draw
    desc: string
        description of what this function is doing; this will be printed by the logger while this function runs

    Returns
    ==================
    lower_xis: list
        list of dimensionless frequency values where a frequency domain waveform enters the detectable region of a detector
    upper_xis: list
        list of dimensionless frequency values where a frequency domain waveform leaves the detectable region of a detector
    '''
    lower_xis = []
    upper_xis = []
    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)
    log.addHandler(ProgressBar())

    for trial in tqdm.tqdm(range(samples), desc = desc):
        roots = []
        injection = prior.sample()
        geometrized_frequency_grid = np.geomspace(0.001,1/np.pi,1000)
        amplitude = np.abs(waveform_generator.frequency_domain_strain(parameters=injection)['plus'])
        M = bilby.gw.conversion.generate_mass_parameters(injection)['total_mass']
        freqs = waveform_generator.frequency_array/float(203025.4467280836/M)
        
        zero_indices = np.where(amplitude==0)[0]
        amplitude = np.delete(amplitude,zero_indices)
        freqs = np.delete(freqs,zero_indices)
        
        effective_amplitude = np.interp(geometrized_frequency_grid,freqs,2*amplitude*np.sqrt(freqs)*np.sqrt(float(203025.4467280836/M)))
        ASD = np.interp(geometrized_frequency_grid,asd_data[:,0]/float(203025.4467280836/M),asd_data[:,1])
        nodes = np.linspace(0,len(geometrized_frequency_grid)-1,100).astype(int)
        parameters = (effective_amplitude-ASD)[nodes]
        spline = scipy.interpolate.CubicSpline(geometrized_frequency_grid[nodes],parameters)
        roots = spline.roots()
        try:
            if roots[0] > 0 and roots[-1] < 1/np.pi:
                if len(roots)>1:
                    lower_xis.append(roots[0])
                    upper_xis.append(roots[-1])
                else:
                    upper_xis.append(roots[0])
        except:
            pass
    return lower_xis, upper_xis



def gaussian(x, N, mu, sigma):
    '''
    Gaussian function.

    Parameters
    ==================
    x: numpy.ndarray
        input array
    N: float
        scaling constant
    mu: float
        mean of the Gaussian
    sigma: float
        standard deviation of the Gaussian

    Returns
    ==================
    output_array: numpy.ndarray
        output Gaussian corresponding to the input array
    '''
    output_array = N*np.exp(-0.5*((x - mu)/sigma)**2)
    
    return output_array



def gaussian_parameters_from_A_ASD_solutions(lower_xis,upper_xis,xi_max):
    '''
    Computes the gaussian fits for the xi_0 and delta_xi_tilde priors. To be used in conjunction with GWCorrect.wfu.utils.A_ASD_solutions

    Parameters
    ==================
    lower_xis: list
        lower_xis from GWCorrect.wfu.utils.A_ASD_solutions
    upper_xis: list
        upper_xis from GWCorrect.wfu.utils.A_ASD_solutions
    xi_max: float
        upper bound on the dimensionless frequency band

    Returns
    ==================
    mu1: float
        mean of the lower_xis distribution
    sigma1: float
        standard deviation of the lower_xis distribution
    mu2: float
        mean of the upper_xis distribution
    sigma2: float
        standard deviation of the upper_xis distribution
    mu: float
        mean of the delta_xi_tilde distribution
    sigma: float
        standard deviation of the delta_xi_tilde distribution
    '''
    lower=lower_xis.copy()
    upper=upper_xis.copy()
    while len(lower) != len(upper):
        upper.pop(0)

    delta_xi_tilde = (np.array(upper)-np.array(lower))/(xi_max-np.array(lower))

    counts, edges = np.histogram(lower_xis, bins=100, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    popt, pcov = scipy.optimize.curve_fit(gaussian, centers, counts, p0=[1, 0.01, 0.01])
    _, mu1, sigma1 = popt

    counts, edges = np.histogram(upper_xis, bins=100, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    popt, pcov = scipy.optimize.curve_fit(gaussian, centers, counts, p0=[1, 0.1, 0.01])
    _, mu2, sigma2 = popt

    counts, edges = np.histogram(delta_xi_tilde, bins=100, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    popt, pcov = scipy.optimize.curve_fit(gaussian, centers, counts, p0=[1, 0.3, 0.1])
    _, mu, sigma = popt
    
    return mu1, np.abs(sigma1), mu2, np.abs(sigma2), mu, np.abs(sigma)
