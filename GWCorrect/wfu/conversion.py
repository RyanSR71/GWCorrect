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



def dimensionless_frequency_nodes_from_delta_xis(xi_0,xi_n,delta_xis):
    '''
    Generates all frequency nodes in dimensionless frequency

    Parameters
    ==================
    xi_0: float
        position of the zeroth frequency node in dimensionless frequency
    xi_n: float
        position of the nth (last) frequency node in dimensionless frequency
    delta_xis: list or numpy.ndarray
        list of the inner frequency node parameters; there must be at least one

    Returns
    ==================
    dimensionless_frequency_nodes: numpy.ndarray
        array of dimensionless_frequency_nodes
    '''
    n = len(delta_xis)+1
    if n == 2:
        xi_1 = xi_0**(1-delta_xis[0])*xi_n**(delta_xis[0])
        dimensionless_frequency_nodes = np.array([xi_0,xi_1,xi_n])
    else:
        xi_matrix = np.zeros((n-1,n-1))
        for l in range(0,n-1):
            for k in range(0,n-1):
                if k == l:
                    xi_matrix[k,l] = 1
                elif k == l-1:
                    xi_matrix[k,l] = -delta_xis[k]
                elif k == l+1:
                    xi_matrix[k,l] = delta_xis[k]-1
        column = np.zeros(n-1)
        column[0] = np.log(xi_0)*(1-delta_xis[0])
        column[-1] = np.log(xi_n)*delta_xis[-1]
        xis = np.exp(np.linalg.solve(xi_matrix,column))
        dimensionless_frequency_nodes = np.array([xi_0]+list(xis)+[xi_n])
        
    return dimensionless_frequency_nodes



def conversion(input_parameters,**kwargs):
    '''
    Conversion function to generate constrained parameters: total_mass, epsilon_f, and epsilon_alpha

    Parameters
    ==================
    input_parameters: dict
        dictionary of BBH source and correction parameters
    n: int, optional
        number of frequency nodes excluding zero; required for epsilon_f and epsilon_alpha calculations
        default: None
    fixed_total_mass: float, optional
        total_mass value for epsilon_f and epsilon_alpha calculations if the total mass is not being sampled in the prior
        default: None
    sigma_A_spline: scipy.interpolate._cubic.CubicSpline
        scipy cubic spline object encoding the standard deviation of the amplitude correction prior; required for epsilon_alpha calculation
        default: None
    gamma: float
        smoothing parameter for the smooth interpolation function
        default: 0.025
    duration: float
        observation time duration of the gravitational-wave signal; required for epsilon_f calculation
        default: None

    Returns
    ==================
    parameters: dict
        dictionary of input parameters plus total mass, epsilon_f, and/or epsilon_alpha; this depends on kwargs
    '''
    n = kwargs.get('n',None)
    fixed_total_mass = kwargs.get('total_mass',None)
    sigma_A_spline = kwargs.get('sigma_A_spline',None)
    gamma = kwargs.get('gamma',0.025)
    duration = kwargs.get('duration',None)
    
    parameters = input_parameters.copy()
    keys = list(parameters.keys())
    if fixed_total_mass is None:
        parameters['total_mass'] = bilby.gw.conversion.generate_mass_parameters(parameters)['total_mass']
    else:
        if hasattr(parameters[keys[0]],'__len__'):
            parameters['total_mass'] = np.array([fixed_total_mass]*len(parameters[keys[0]]))
        else:
            parameters['total_mass'] = fixed_total_mass
    if n is not None:
        if hasattr(parameters[keys[0]],'__len__'):
            epsilon_f = []
            epsilon_alpha = []
            for i in range(len(parameters[keys[0]])):
                delta_f = 0
                delta_xis = [parameters[f'delta_xi_{k}'][i] for k in range(1,n)]
                total_mass = parameters['total_mass'][i]
                xi_0 = parameters['xi_0'][i]
                xi_n = parameters['xi_n'][i]
                dimensionless_frequency_nodes = dimensionless_frequency_nodes_from_delta_xis(xi_0,xi_n,delta_xis)
                delta_f_0n = wfu.conversion.convert_to_SI_frequency(xi_n-xi_0,total_mass)
                if delta_f_0n < 4*n/duration:
                    delta_f += 1
                for k in range(1,n+1):
                    delta_f_jk = wfu.conversion.convert_to_SI_frequency(dimensionless_frequency_nodes[k]-dimensionless_frequency_nodes[k-1],total_mass)
                    if delta_f_jk < 4/duration:
                        delta_f += 1
                epsilon_f.append(delta_f)
                
                if sigma_A_spline is not None:
                    alpha_tildes = [0]+[parameters[f'alpha_tilde_{k}'][i] for k in range(1,n+1)]
                    alphas = np.array(alpha_tildes)*sigma_A_spline(dimensionless_frequency_nodes)
                    if any(alphas[k] < -1 for k in range(n+1)):
                        epsilon_alpha.append(1)
                    else:
                        epsilon_alpha.append(0)
                else:
                    epsilon_alpha.append(0)
                    
            epsilon_f = np.array(epsilon_f)
            epsilon_alpha = np.array(epsilon_alpha)
            
        else:
            epsilon_f = 0
            delta_xis = [parameters[f'delta_xi_{k}'] for k in range(1,n)]
            total_mass = parameters['total_mass']
            xi_0 = parameters['xi_0']
            xi_n = parameters['xi_n']
            dimensionless_frequency_nodes = dimensionless_frequency_nodes_from_delta_xis(xi_0,xi_n,delta_xis)
            delta_f_0n = wfu.conversion.convert_to_SI_frequency(xi_n-xi_0,total_mass)
            if delta_f_0n < 4*n/duration:
                epsilon_f += 1
            for k in range(1,n+1):
                delta_f_jk = wfu.conversion.convert_to_SI_frequency(dimensionless_frequency_nodes[k]-dimensionless_frequency_nodes[k-1],total_mass)
                if delta_f_jk < 4/duration:
                    epsilon_f += 1
            
            epsilon_alpha = 0
            if sigma_A_spline is not None:
                alpha_tildes = [0]+[parameters[f'alpha_tilde_{k}'] for k in range(1,n+1)]
                alphas = np.array(alpha_tildes)*sigma_A_spline(dimensionless_frequency_nodes)
                if any(alphas[k] < -1 for k in range(n+1)):
                    epsilon_alpha += 1
                    
        parameters['epsilon_f'] = epsilon_f
        parameters['epsilon_alpha'] = epsilon_alpha
        
    return parameters
