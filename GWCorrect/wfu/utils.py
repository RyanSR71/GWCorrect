import numpy as np
import bilby
import scipy
import tqdm
import logging
from .conversion import convert_to_dimensionless_frequency, convert_to_SI_frequency



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



def amplitude_difference_from_parameters(frequency_array,dimensionless_frequency_nodes,alpha_tildes,sigma_A_spline,mass_1,mass_2,gamma):
    '''
    Generates the amplitude correction from the frequency nodes and amplitude correction parameters

    Parameters
    ==================
    frequency_array: numpy.ndarray
        array of SI frequency points over which to generate the amplitude correction
    dimensionless_frequency_nodes: numpy.ndarray
        array of frequency nodes in dimensionless frequency nodes
    alpha_tildes: list
        list of rescaled amplitude correction parameters
    sigma_A_spline: scipy.interpolate._cubic.CubicSpline
        scipy cubic spline object encoding the total standard deviation of the amplitude correction prior
    mass_1: float
        primary mass in solar masses
    mass_2: float
        secondary mass in solar masses
    gamma: float
        smoothing parameter for the smooth interpolation function
    '''
    n = len(alpha_tildes)-1
    total_mass = mass_1+mass_2
    dimensionless_frequency_array = np.linspace(dimensionless_frequency_nodes[0],dimensionless_frequency_nodes[-1],len(frequency_array))
    sigma_A = sigma_A_spline(dimensionless_frequency_nodes)

    amplitude_difference_in_dimensionless_frequency = smooth_interpolation(dimensionless_frequency_array,dimensionless_frequency_nodes,alpha_tildes*sigma_A,gamma)
    amplitude_difference = np.interp(frequency_array,convert_to_SI_frequency(dimensionless_frequency_array,total_mass),amplitude_difference_in_dimensionless_frequency)

    return amplitude_difference



def phase_difference_from_parameters(frequency_array,dimensionless_frequency_nodes,phi_tildes,sigma_phi_spline,mass_1,mass_2,gamma):
    '''
    Generates the phase correction from the frequency nodes and phase correction parameters

    Parameters
    ==================
    frequency_array: numpy.ndarray
        array of SI frequency points over which to generate the phase correction
    dimensionless_frequency_nodes: numpy.ndarray
        array of frequency nodes in dimensionless frequency nodes
    phi_tildes: list
        list of rescaled phase correction parameters
    sigma_phi_spline: scipy.interpolate._cubic.CubicSpline
        scipy cubic spline object encoding the total standard deviation of the phase correction prior
    mass_1: float
        primary mass in solar masses
    mass_2: float
        secondary mass in solar masses
    gamma: float
        smoothing parameter for the smooth interpolation function
    '''
    n = len(phi_tildes)-1
    total_mass = mass_1+mass_2
    dimensionless_frequency_array = np.linspace(dimensionless_frequency_nodes[0],dimensionless_frequency_nodes[-1],len(frequency_array))
    sigma_phi = sigma_phi_spline(dimensionless_frequency_nodes)

    phase_difference_in_dimensionless_frequency = smooth_interpolation(dimensionless_frequency_array,dimensionless_frequency_nodes,phi_tildes*sigma_phi,gamma)
    phase_difference = np.interp(frequency_array,convert_to_SI_frequency(dimensionless_frequency_array,total_mass),phase_difference_in_dimensionless_frequency)

    return phase_difference



def A_ASD_solutions(waveform_generator,asd_data,prior,minimum_frequency,maximum_frequency,samples,desc):
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
    minimum_frequency: float
        lower bound on the SI frequency band
    maximum_frequency: float
        upper bound on the SI frequency band
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
        freqs = convert_to_dimensionless_frequency(waveform_generator.frequency_array,M)
        
        zero_indices = np.where(amplitude==0)[0]
        amplitude = np.delete(amplitude,zero_indices)
        freqs = np.delete(freqs,zero_indices)
        
        effective_amplitude = np.interp(geometrized_frequency_grid,freqs,2*amplitude*np.sqrt(waveform_generator.frequency_array))
        ASD = np.interp(geometrized_frequency_grid,convert_to_dimensionless_frequency(asd_data[:,0],M),asd_data[:,1])
        nodes = np.linspace(0,len(geometrized_frequency_grid)-1,100).astype(int)
        parameters = (effective_amplitude-ASD)[nodes]
        spline = scipy.interpolate.CubicSpline(geometrized_frequency_grid[nodes],parameters)
        roots = spline.roots()
        try:
            if roots[0] > 0 and roots[-1] < 1/np.pi:
                if len(roots)>1:
                    for root in roots:
                        if root > convert_to_dimensionless_frequency(minimum_frequency,M):
                            lower_xis.append(root)
                            break
                    for root in roots[::-1]:
                        if root < convert_to_dimensionless_frequency(maximum_frequency,M):
                            upper_xis.append(root)
                            break
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



def gaussian_parameters_from_A_ASD_solutions(lower_xis,upper_xis):
    '''
    Computes the gaussian fits for the xi_0 and delta_xi_tilde priors. To be used in conjunction with GWCorrect.wfu.utils.A_ASD_solutions

    Parameters
    ==================
    lower_xis: list
        lower_xis from GWCorrect.wfu.utils.A_ASD_solutions
    upper_xis: list
        upper_xis from GWCorrect.wfu.utils.A_ASD_solutions

    Returns
    ==================
    mu_0: float
        mean of the lower_xis distribution
    sigma_0: float
        standard deviation of the lower_xis distribution
    mu_n: float
        mean of the upper_xis distribution
    sigma_n: float
        standard deviation of the upper_xis distribution
    '''

    counts, edges = np.histogram(lower_xis, bins=100, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    popt, pcov = scipy.optimize.curve_fit(gaussian, centers, counts, p0=[1, 0.01, 0.01])
    _, mu_0, sigma_0 = popt

    counts, edges = np.histogram(upper_xis, bins=100, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    popt, pcov = scipy.optimize.curve_fit(gaussian, centers, counts, p0=[1, 0.1, 0.01])
    _, mu_n, sigma_n = popt
    
    return mu_0, np.abs(sigma_0), mu_n, np.abs(sigma_n)
