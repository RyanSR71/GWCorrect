def xi_0_upper_bound(n, **kwargs):
    xi_min = kwargs.get('xi_min',0.018)
    xi_max = kwargs.get('xi_max',1/np.pi)
    def f(x):
        return x**(1-n) + (xi_min**(1-n)/(xi_max-xi_min))*x - ((xi_max*xi_min**(1-n))/(xi_max-xi_min))
    # we need to have the lower bound slightly larger than xi_min, so we multiply by 1.000001 arbitrarily
    root = scipy.optimize.brentq(f, 1.000001*xi_min, xi_max)
    return root



def delta_xi_tilde_lower_bound(n,minimum_frequency,duration,**kwargs):
    xi_min = kwargs.get('xi_min',0.018)
    xi_max = kwargs.get('xi_max',1/np.pi)
    
    minimum = (xi_min/(xi_max-xi_min))*(4/(duration*minimum_frequency))**n
    
    return minimum



def td_waveform(fd_waveform,sampling_frequency):
    reversed_fd_waveform = fd_waveform[::-1]
    total_fd_strain = np.concatenate((fd_waveform,np.conjugate(reversed_fd_waveform[1:-1])))
    td_waveform = sampling_frequency*np.real(np.fft.ifft(total_fd_strain))
    return td_waveform



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



def full_conversion(parameters,minimum_frequency,maximum_frequency,xi_max,n,sigma_dphi_spline):
    '''
    Conversion function to generate the total mass, time shift, phase shift, and frequency node constraints from a set of parameters;
    Necessary for the use of the constraint priors
    
    Parameters
    ==================
    parameters: dict
        dictionary of binary black hole parameters
    minimum_frequency: float
        lower bound on the frequency band (Hz)
    maximum_frequency: float
        upper bound on the frequency band (Hz)
    xi_max: float
        absolute upper bound on dimensionless frequency
    n: int
        number of frequency nodes (excluding the zeroth)
    sigma_dphi_spline: scipy.interpolate.CubicSpline
        spline object representing the standard deviation of phase model differences in dimensionless frequency
    
    Returns
    ==================
    parameters: dict
        input parameters, but with the constraint parameters added
    '''
    total_mass = bilby.gw.conversion.generate_mass_parameters(parameters)['total_mass']
    
    if hasattr(total_mass,"__len__"):
        row_number = len(total_mass)
    else:
        row_number = 1
    
    constraint_matrix = np.zeros((row_number,2*n+12))

    constraint_matrix[:,0] = total_mass   
    constraint_matrix[:,1] = parameters['xi_0']
    constraint_matrix[:,2] = parameters['delta_xi_tilde']
    for k in range(1,n+1):
        constraint_matrix[:,n+7+k] = parameters[f'dphi_{k}']

    for i in range(row_number):
        constraint_matrix[i,3] = minimum_frequency
        constraint_matrix[i,n+5] = maximum_frequency
        for j in range(4,n+5):
            constraint_matrix[i,j] = (203025.4467280836*constraint_matrix[i,1]/constraint_matrix[i,0])*(1+((xi_max-constraint_matrix[i,1])/constraint_matrix[i,1])*constraint_matrix[i,2])**((j-4)/n)
            constraint_matrix[i,j+n+3] *= sigma_dphi_spline(4.925490947641267e-06*constraint_matrix[i,0]*constraint_matrix[i,j])
        constraint_matrix[i,2*n+8] = constraint_matrix[i,2*n+7]
        frequency_nodes = constraint_matrix[i][3:n+6]
        dphis = constraint_matrix[i][n+6:2*n+9]
        constraint_matrix[i][2*n+9],constraint_matrix[i][2*n+10] = np.polyfit(frequency_nodes,dphis,1)
        constraint_matrix[i][2*n+9] *= 1/(2*np.pi)
        constraint_matrix[i][2*n+11] = 1/(constraint_matrix[i,5]-constraint_matrix[i,4])

    t_dphi = constraint_matrix[:,2*n+9]
    phi_dphi = constraint_matrix[:,2*n+10]
    delta_t_01 = constraint_matrix[:,2*n+11]
    
    if row_number == 1:
        parameters['total_mass'] = float(total_mass)
        parameters['t_dphi'] = float(t_dphi)
        parameters['phi_dphi'] = float(phi_dphi)
        parameters['delta_t_01'] = float(delta_t_01)
    else:    
        parameters['total_mass'] = total_mass
        parameters['t_dphi'] = t_dphi
        parameters['phi_dphi'] = phi_dphi
        parameters['delta_t_01'] = delta_t_01
    
    return parameters



class BasicCorrectionModel(object):
    '''
    Modified WaveformGenerator object from bilby.gw to include waveform uncertainty corrections in the strain calculation
    To sample waveform uncertainty, include all relevant "alpha" and "phi" parameters in the prior.
    Note: make sure the number of alphas, phis, and waveform_uncertainty_nodes are the same
    
    New Parameters
    ==================
   frequency_nodes: numpy.ndarray, optional
        array of frequency nodes to be used in generating the dA and dphi splines
        default: None
    correct_amplitude: bool, optional
        if True, the waveform generator will attempt to pull dA parameters from the parameter dictionary (either an injection or the prior)
        default: None
    correct_phase: bool, optional
        if True, the waveform generator will attempt to pull dphi parameters from the parameter dictionary (either an injection or the prior)
        default: None
    dimensionless: bool, optional
        whether or not to perfom the waveform correction in dimensionless frequency or not
        default: True
    '''
    duration = PropertyAccessor('_times_and_frequencies', 'duration')
    sampling_frequency = PropertyAccessor('_times_and_frequencies', 'sampling_frequency')
    start_time = PropertyAccessor('_times_and_frequencies', 'start_time')
    frequency_array = PropertyAccessor('_times_and_frequencies', 'frequency_array')
    time_array = PropertyAccessor('_times_and_frequencies', 'time_array')
    def __init__(self, duration=None, sampling_frequency=None, start_time=0, frequency_domain_source_model=None,
                 time_domain_source_model=None, parameters=None,
                 frequency_nodes=None,correct_amplitude=False,correct_phase=False,
                 dimensionless=True,parameter_conversion=None,
                 waveform_arguments=None):
        self._times_and_frequencies = CoupledTimeAndFrequencySeries(duration=duration,
                                                                    sampling_frequency=sampling_frequency,
                                                                    start_time=start_time)
        self.frequency_domain_source_model = frequency_domain_source_model
        self.time_domain_source_model = time_domain_source_model
        self.source_parameter_keys = self.__parameters_from_source_model()
        self.dimensionless=dimensionless
        self.correct_amplitude=correct_amplitude
        self.correct_phase=correct_phase
        
        if parameter_conversion is None:
            self.parameter_conversion = convert_to_lal_binary_black_hole_parameters
        else:
            self.parameter_conversion = parameter_conversion
        if waveform_arguments is not None:
            self.waveform_arguments = waveform_arguments
        else:
            self.waveform_arguments = dict()
        if frequency_nodes is not None:
            self.frequency_nodes = frequency_nodes
        else:
            self.frequency_nodes = None
        if isinstance(parameters, dict):
            self.parameters = parameters
        self._cache = dict(parameters=None, waveform=None, model=None)
        
        utils.logger.info(
            "Waveform generator initiated with\n"
            "  frequency_domain_source_model: {}\n"
            "  time_domain_source_model: {}\n"
            "  parameter_conversion: {}"
            .format(utils.get_function_path(self.frequency_domain_source_model),
                    utils.get_function_path(self.time_domain_source_model),
                    utils.get_function_path(self.parameter_conversion))
        )
        
    def __repr__(self):
        if self.frequency_domain_source_model is not None:
            fdsm_name = self.frequency_domain_source_model.__name__
        else:
            fdsm_name = None
        if self.time_domain_source_model is not None:
            tdsm_name = self.time_domain_source_model.__name__
        else:
            tdsm_name = None
        if self.parameter_conversion is None:
            param_conv_name = None
        else:
            param_conv_name = self.parameter_conversion.__name__
        return self.__class__.__name__ + '(duration={}, sampling_frequency={}, start_time={}, ' \
                                         'frequency_domain_source_model={}, time_domain_source_model={}, ' \
                                         'parameter_conversion={}, ' \
                                         'frequency_nodes={}, ' \
                                         'dimensionless={}, ' \
                                         'correct_amplitude={}, ' \
                                         'correct_phase={}, ' \
                                         'waveform_arguments={})'\
            .format(self.duration, self.sampling_frequency, self.start_time, fdsm_name, tdsm_name,
                    param_conv_name, self.frequency_nodes, self.dimensionless, self.correct_amplitude, self.correct_phase, self.waveform_arguments)
    
    def frequency_domain_strain(self, parameters=None):
        return self._calculate_strain(model=self.frequency_domain_source_model,
                                      model_data_points=self.frequency_array,
                                      parameters=parameters,
                                      transformation_function=utils.nfft,
                                      transformed_model=self.time_domain_source_model,
                                      transformed_model_data_points=self.time_array,
                                      frequency_nodes=self.frequency_nodes,
                                      dimensionless=self.dimensionless,
                                      correct_amplitude=self.correct_amplitude,
                                      correct_phase=self.correct_phase,
                                      )

    def time_domain_strain(self,parameters=None):
        fd_model_strain = self._calculate_strain(model=self.frequency_domain_source_model,
                                      model_data_points=self.frequency_array,
                                      parameters=parameters,
                                      transformation_function=utils.nfft,
                                      transformed_model=self.time_domain_source_model,
                                      transformed_model_data_points=self.time_array,
                                      frequency_nodes=self.frequency_nodes,
                                      dimensionless=self.dimensionless,
                                      correct_amplitude=self.correct_amplitude,
                                      correct_phase=self.correct_phase,
                                      )

        model_strain = dict()
        model_strain['plus'] = td_waveform(fd_model_strain['plus'],self.sampling_frequency)
        model_strain['cross'] = td_waveform(fd_model_strain['cross'],self.sampling_frequency)

        return model_strain
    
    def _calculate_strain(self, model, model_data_points, transformation_function, transformed_model,
                          transformed_model_data_points, parameters, frequency_nodes, dimensionless, correct_amplitude, correct_phase):
        if parameters is not None:
            self.parameters = parameters
        if self.parameters == self._cache['parameters'] and self._cache['model'] == model and \
                self._cache['transformed_model'] == transformed_model:
            return self._cache['waveform']
        if model is not None:
            model_strain = self._strain_from_model(model_data_points, model)
        elif transformed_model is not None:
            model_strain = self._strain_from_transformed_model(transformed_model_data_points, transformed_model,
                                                               transformation_function)
        else:
            raise RuntimeError("No source model given")
        
        '''
        The following block performs the waveform uncertainty correction:
        '''
        
        if self.frequency_nodes is not None:
            
            indexes = np.arange(0,len(self.frequency_nodes),1)
            
            if self.dimensionless is True:
                M = bilby.gw.conversion.generate_mass_parameters(parameters)['total_mass']
                frequency_nodes = np.array(self.frequency_nodes)*float(203025.4467280836/M)
            else:
                frequency_nodes = self.frequency_nodes

            try:
                gamma = parameters['gamma']
            except:
                gamma = 0.025
            
            if self.correct_amplitude is True:
                try:
                    alphas = [parameters[f'dA_{i}'] for i in indexes]
                    dA = smooth_interpolation(self.frequency_array,frequency_nodes,alphas,gamma)
                except:
                    raise Exception('Amplitude Correction Failed!')
            else:
                dA = 0

            if self.correct_phase is True:
                try:
                    phis = [parameters[f'dphi_{i}'] for i in indexes]
                    dphi = smooth_interpolation(self.frequency_array,frequency_nodes,phis,gamma)
                except:
                    raise Exception('Phase Correction Failed!')
            else:
                dphi = 0
                    
        else:
            dA = 0
            dphi = 0
            
        model_strain['plus'] = model_strain['plus']*(1+dA)*np.exp(dphi*1j)
        model_strain['cross'] = model_strain['cross']*(1+dA)*np.exp(dphi*1j)
        
        self._cache['waveform'] = model_strain
        self._cache['parameters'] = self.parameters.copy()
        self._cache['model'] = model
        self._cache['transformed_model'] = transformed_model
        self._cache['amplitude_difference'] = dA
        self._cache['phase_difference'] = dphi
                              
        return model_strain
    
    def _strain_from_model(self, model_data_points, model):
        return model(model_data_points, **self.parameters)
    
    def _strain_from_transformed_model(self, transformed_model_data_points, transformed_model, transformation_function):
        transformed_model_strain = self._strain_from_model(transformed_model_data_points, transformed_model)
        if isinstance(transformed_model_strain, np.ndarray):
            return transformation_function(transformed_model_strain, self.sampling_frequency)
        model_strain = dict()
        for key in transformed_model_strain:
            if transformation_function == bilby.utils.nfft:
                model_strain[key], _ = \
                    transformation_function(transformed_model_strain[key], self.sampling_frequency)
            else:
                model_strain[key] = transformation_function(transformed_model_strain[key], self.sampling_frequency)
        return model_strain
    
    @property
    def parameters(self):
        return self.__parameters
    
    @parameters.setter 
    def parameters(self, parameters):
        if not isinstance(parameters, dict):
            raise TypeError('"parameters" must be a dictionary.')
        new_parameters = parameters.copy()
        new_parameters, _ = self.parameter_conversion(new_parameters)
        for key in self.source_parameter_keys.symmetric_difference(new_parameters):
            # preventing waveform uncertainty parameters from being removed
            if self.frequency_nodes is not None:
                indexes = np.arange(0,len(self.frequency_nodes),1)
                if key not in [f'dA_{i}' for i in indexes]+[f'dphi_{i}' for i in indexes]+['gamma_dphi','gamma_dA']:  
                    new_parameters.pop(key)
            else:
                new_parameters.pop(key)
        self.__parameters = new_parameters
        self.__parameters.update(self.waveform_arguments)
        
    def __parameters_from_source_model(self):
        if self.frequency_domain_source_model is not None:
            model = self.frequency_domain_source_model
        elif self.time_domain_source_model is not None:
            model = self.time_domain_source_model
        else:
            raise AttributeError('Either time or frequency domain source '
                                 'model must be provided.')
        return set(utils.infer_parameters_from_function(model))



class GeneralCorrectionModel(object):
    '''
    Modified WaveformGenerator object from bilby.gw to include waveform uncertainty corrections in the strain calculation
    
    New Parameters
    ==================
    correction_arguments: dict, optional
        dictionary containing arguments for the waveform correction
        default: None

        contents:
            correct_amplitude: bool
                whether or not to attempt an amplitude correction
                include all dA parameters in the prior if True
            sigma_dA: numpy.ndarry
                if correct_amplitude is True, this is the array of standard deviations for the dA priors
            correct_phase: bool
                whether or not to attempt a phase correction
                include all dphi parameters in the prior if True
            sigma_dphi: numpy.ndarry
                if correct_phase is True, this is the array of standard deviations for the dphi priors
            nodes: int
                number of frequency nodes desired
            xi_high: float, optional
                absolute lower bound on dimensionless frequency, xi
                default: 1/pi
    '''
    duration = PropertyAccessor('_times_and_frequencies', 'duration')
    sampling_frequency = PropertyAccessor('_times_and_frequencies', 'sampling_frequency')
    start_time = PropertyAccessor('_times_and_frequencies', 'start_time')
    frequency_array = PropertyAccessor('_times_and_frequencies', 'frequency_array')
    time_array = PropertyAccessor('_times_and_frequencies', 'time_array')
    def __init__(self, duration=None, sampling_frequency=None, start_time=0, frequency_domain_source_model=None,
                 time_domain_source_model=None, parameters=None,
                 correction_arguments=None,parameter_conversion=None,
                 waveform_arguments=None):
        self._times_and_frequencies = CoupledTimeAndFrequencySeries(duration=duration,
                                                                    sampling_frequency=sampling_frequency,
                                                                    start_time=start_time)
        self.frequency_domain_source_model = frequency_domain_source_model
        self.time_domain_source_model = time_domain_source_model
        self.source_parameter_keys = self.__parameters_from_source_model()
        self.correction_arguments=correction_arguments
        
        if parameter_conversion is None:
            self.parameter_conversion = convert_to_lal_binary_black_hole_parameters
        else:
            self.parameter_conversion = parameter_conversion
        if waveform_arguments is not None:
            self.waveform_arguments = waveform_arguments
        else:
            self.waveform_arguments = dict()
        if self.correction_arguments is not None:
            self.correction_arguments = correction_arguments
        else:
            self.correction_arguments = None
        if isinstance(parameters, dict):
            self.parameters = parameters
        self._cache = dict(parameters=None, waveform=None, model=None)
        
        utils.logger.info(
            "Waveform generator initiated with\n"
            "  frequency_domain_source_model: {}\n"
            "  time_domain_source_model: {}\n"
            "  parameter_conversion: {}"
            .format(utils.get_function_path(self.frequency_domain_source_model),
                    utils.get_function_path(self.time_domain_source_model),
                    utils.get_function_path(self.parameter_conversion))
        )
        
    def __repr__(self):
        if self.frequency_domain_source_model is not None:
            fdsm_name = self.frequency_domain_source_model.__name__
        else:
            fdsm_name = None
        if self.time_domain_source_model is not None:
            tdsm_name = self.time_domain_source_model.__name__
        else:
            tdsm_name = None
        if self.parameter_conversion is None:
            param_conv_name = None
        else:
            param_conv_name = self.parameter_conversion.__name__
        return self.__class__.__name__ + '(duration={}, sampling_frequency={}, start_time={}, ' \
                                         'frequency_domain_source_model={}, time_domain_source_model={}, ' \
                                         'parameter_conversion={}, ' \
                                         'correction_arguments={}, ' \
                                         'waveform_arguments={})'\
            .format(self.duration, self.sampling_frequency, self.start_time, fdsm_name, tdsm_name,
                    param_conv_name, self.correction_arguments, self.waveform_arguments)
    
    def frequency_domain_strain(self, parameters=None):
        return self._calculate_strain(model=self.frequency_domain_source_model,
                                      model_data_points=self.frequency_array,
                                      parameters=parameters,
                                      transformation_function=utils.nfft,
                                      transformed_model=self.time_domain_source_model,
                                      transformed_model_data_points=self.time_array,
                                      correction_arguments=self.correction_arguments,
                                      )

    def time_domain_strain(self,parameters=None):
        fd_model_strain = self._calculate_strain(model=self.frequency_domain_source_model,
                                      model_data_points=self.frequency_array,
                                      parameters=parameters,
                                      transformation_function=utils.nfft,
                                      transformed_model=self.time_domain_source_model,
                                      transformed_model_data_points=self.time_array,
                                      correction_arguments=self.correction_arguments,
                                      )

        model_strain = dict()
        model_strain['plus'] = td_waveform(fd_model_strain['plus'],self.sampling_frequency)
        model_strain['cross'] = td_waveform(fd_model_strain['cross'],self.sampling_frequency)

        return model_strain
    
    def _calculate_strain(self, model, model_data_points, transformation_function, transformed_model,
                          transformed_model_data_points, parameters, correction_arguments):
        if parameters is not None:
            self.parameters = parameters
        if self.parameters == self._cache['parameters'] and self._cache['model'] == model and \
                self._cache['transformed_model'] == transformed_model:
            return self._cache['waveform']
        if model is not None:
            model_strain = self._strain_from_model(model_data_points, model)
        elif transformed_model is not None:
            model_strain = self._strain_from_transformed_model(transformed_model_data_points, transformed_model,
                                                               transformation_function)
        else:
            raise RuntimeError("No source model given")
        
        '''
        The following block performs the waveform uncertainty correction:
        '''
                              
        try:
            gamma = parameters['gamma']
        except:
            gamma = 0.025

        try:
            xi_high = correction_arguments['xi_high']
        except:
            xi_high = 1/np.pi

        M = bilby.gw.conversion.generate_mass_parameters(parameters)['total_mass']
        indexes = np.arange(0,correction_arguments['nodes']+1,1)
                              
        if correction_arguments['correct_amplitude'] is True:
            sigma_dA = correction_arguments['sigma_dA']
            xi_up = parameters['xi_0']*(1+((xi_high-parameters['xi_0'])/(parameters['xi_0']))*parameters['delta_xi_tilde'])
            dA_frequency_nodes,dA_coeffs = variable_prior(sigma_dA,correction_arguments['nodes'],
                                                      parameters['xi_0'], xi_up)
            dA_frequency_nodes *= float(203025.4467280836/M)
            try:
                prior_alphas = np.array([parameters[f'dA_{i}'] for i in indexes])
                alphas = prior_alphas*dA_coeffs
                dA = smooth_interpolation(self.frequency_array,dA_frequency_nodes,alphas,gamma)
            except:
                raise Exception('Amplitude Correction Failed!')
        else:
            dA = 0
        if correction_arguments['correct_phase'] is True:
            sigma_dphi = correction_arguments['sigma_dphi']
            xi_up = parameters['xi_0']*(1+((xi_high-parameters['xi_0'])/(parameters['xi_0']))*parameters['delta_xi_tilde'])
            dphi_frequency_nodes,dphi_coeffs = variable_prior(sigma_dphi,correction_arguments['nodes'],
                                                            parameters['xi_0'], xi_up)
            dphi_frequency_nodes *= float(203025.4467280836/M)
            try:
                prior_phis = np.array([parameters[f'dphi_{i}'] for i in indexes])
                phis = prior_phis*dphi_coeffs
                dphi = smooth_interpolation(self.frequency_array,dphi_frequency_nodes,phis,gamma)
            except:
                raise Exception('Phase Correction Failed!')
        else:
            dphi = 0

        model_strain['plus'] = model_strain['plus']*(1+dA)*np.exp(dphi*1j)
        model_strain['cross'] = model_strain['cross']*(1+dA)*np.exp(dphi*1j)
        
        self._cache['waveform'] = model_strain
        self._cache['parameters'] = self.parameters.copy()
        self._cache['model'] = model
        self._cache['transformed_model'] = transformed_model
        self._cache['amplitude_difference'] = dA
        self._cache['phase_difference'] = dphi
                              
        return model_strain
    
    def _strain_from_model(self, model_data_points, model):
        return model(model_data_points, **self.parameters)
    
    def _strain_from_transformed_model(self, transformed_model_data_points, transformed_model, transformation_function):
        transformed_model_strain = self._strain_from_model(transformed_model_data_points, transformed_model)
        if isinstance(transformed_model_strain, np.ndarray):
            return transformation_function(transformed_model_strain, self.sampling_frequency)
        model_strain = dict()
        for key in transformed_model_strain:
            if transformation_function == utils.nfft:
                model_strain[key], _ = \
                    transformation_function(transformed_model_strain[key], self.sampling_frequency)
            else:
                model_strain[key] = transformation_function(transformed_model_strain[key], self.sampling_frequency)
        return model_strain
    
    @property
    def parameters(self):
        return self.__parameters
    
    @parameters.setter 
    def parameters(self, parameters):
        if not isinstance(parameters, dict):
            raise TypeError('"parameters" must be a dictionary.')
        new_parameters = parameters.copy()
        new_parameters, _ = self.parameter_conversion(new_parameters)
        for key in self.source_parameter_keys.symmetric_difference(new_parameters):
            # preventing waveform uncertainty parameters from being removed
            indexes = np.arange(0,self.correction_arguments['nodes']+1,1)
            if key not in [f'dA_{i}' for i in indexes]+[f'dphi_{i}' for i in indexes]+['gamma','xi_0','delta_xi_tilde']:  
                new_parameters.pop(key)
        self.__parameters = new_parameters
        self.__parameters.update(self.waveform_arguments)
        
    def __parameters_from_source_model(self):
        if self.frequency_domain_source_model is not None:
            model = self.frequency_domain_source_model
        elif self.time_domain_source_model is not None:
            model = self.time_domain_source_model
        else:
            raise AttributeError('Either time or frequency domain source '
                                 'model must be provided.')
        return set(utils.infer_parameters_from_function(model))
