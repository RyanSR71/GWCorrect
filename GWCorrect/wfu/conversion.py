import numpy as np
import scipy
import bilby

def convert_to_dimensionless_frequency(frequencies,total_mass):
    '''
    Converts SI frequency to mass-weighted dimensionless frequency.

    Parameters
    ==================
    frequencies: float or numpy.ndarray
        SI frequency values
    total_mass: float
        total binary mass to weight the frequencies

    Returns
    ==================
    dimensionless_frequency: float or numpy.ndarray
        dimensionless frequencies of the same shape as the input frequencies
    '''
    return 0.000004925490947641267*total_mass*frequencies

def convert_to_SI_frequency(dimensionless_frequencies,total_mass):
    '''
    Converts mass-weighted dimensionless frequency to SI frequency.

    Parameters
    ==================
    dimensionless_frequencies: float or numpy.ndarray
        mass-weighted dimensionless frequency values
    total_mass: float
        total binary mass to weight the frequencies

    Returns
    ==================
    SI_frequency: float or numpy.ndarray
        SI frequencies of the same shape as the input dimensionless frequencies
    '''
    return 203025.4467280836*dimensionless_frequencies/total_mass

def dimensionless_frequency_nodes_from_xi_0_delta_xi_tilde(xi_0,delta_xi_tilde,n,xi_max):
    '''
    Converts xi_0 and delta_xi_tilde to frequency nodes in mass-weighted dimensionless frequency.

    Parameters
    ==================
    xi_0: float
        starting frequency node parameter
    delta_xi_tilde: float
        frequency_node_spacing_parameter
    n: int
        number of frequency nodes excluding 0
    xi_max: float
        upper bound on the dimensionless frequency band
    '''
    dimensionless_frequency_nodes = [xi_0*(1+((xi_max-xi_0)/xi_0)*delta_xi_tilde)**(k/n) for k in range(0,n+1)]
    return dimensionless_frequency_nodes

def delta_t_01_from_frequency_node_parameters(xi_0,delta_xi_tilde,total_mass,n,xi_max):
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


def delta_t_01_posterior(result,n,xi_max):
    '''
    Computes delta_t_01 from a bilby result object. Used to add delta_t_01 to the posterior.

    Parameters
    ==================
    result: bilby.core.result.Result
        bilby result object
    n: int
        number of frequency nodes excluding 0
    xi_max: float
        upper bound on the dimensionless frequency band

    Returns
    ==================
    delta_t_01_posterior: numpy.ndarray
        delta_t_01 posterior samples
    '''
    delta_t_01_posterior = delta_t_01_from_frequency_node_parameters(result.posterior['xi_0'],result.posterior['delta_xi_tilde'],result.posterior['total_mass'],n,xi_max)
    return delta_t_01_posterior



def epsilon_alpha_conversion(dAs,xi_0,delta_xi_tilde,xi_max,sigma_dA_spline):
    '''
    Computes epsilon_alpha up to 1.

    Parameters
    ==================
    dAs: list or numpy.ndarray
        list of dA parameters including dA_0
    xi_0: float
        starting frequency node parameter
    delta_xi_tilde: float
        frequency node spacing parameter
    xi_max: float
        upper bound on the dimensionless frequency band
    sigma_dA_spline: scipy.interpolate._cubic.CubicSpline
        scipy cubic spline object encoding the standard deviation of the dA prior

    Returns
    ==================
    epsilon_alpha: float
        amplitude error parameter
    '''
    n = len(dAs)-1
    dimensionless_frequency_nodes = [xi_0*(1+((xi_max-xi_0)/(xi_0))*delta_xi_tilde)**k/n for k in range(0,n+1)]
    true_dAs = [dAs[k]*sigma_dA_spline(dimensionless_frequency_nodes[k]) for k in range(0,n+1)] 
    counter = n
    if any(x < 0 for x in list(np.array(true_dAs)+1)):
        counter -= 1
    epsilon_alpha = n-counter
    
    return epsilon_alpha
