import numpy as np
import bilby
import random
import time
import sys
import scipy
import lal
import tqdm
import logging
from bilby.core import utils
from bilby.core.series import CoupledTimeAndFrequencySeries
from bilby.core.utils import PropertyAccessor
from bilby.gw.conversion import convert_to_lal_binary_neutron_star_parameters



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
    new_data = data.copy()
    
    r = int(gamma*len(data))
    if r != 0:
        lower_index = int(gamma*len(data))
        upper_index = int((1-gamma)*len(data))
        for i in range(lower_index,upper_index):
            new_data[i] = (1/(2*r))*np.sum(data[i-r:i+r])
    
    output = np.interp(full_grid,temp_grid,new_data)
    if np.abs(output[0]) > 0:
        output -= output[0]

    return output



def GC_waveform_correction(frequency_array,xi_0,delta_xi_tilde,dAs,dphis,sigma_dA_spline,sigma_dphi_spline,mass_1,mass_2,xi_max,gamma):
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



def delta_t_01(xi_0,delta_xi_tilde,total_mass,n,xi_max):
    '''
    Computes delta_t_01.

    Parameters
    ==================
    xi_0: float
        starting frequency node parameter
    delta_xi_tilde: float
        frequency node spacing parameter
    total_mass: float
        total system mass (m_1+m_2)
    n: int
        number of frequency nodes excluding 0
    xi_max: float
        upper bound on the dimensionless frequency band

    Returns
    ==================
    delta_t_01: float
        delta_t_01 value
    '''
    delta_t_01 = 0.000004925490947641267*n*total_mass/((xi_max-xi_0)*delta_xi_tilde)

    return delta_t_01



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
