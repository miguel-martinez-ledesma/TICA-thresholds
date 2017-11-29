%% 
% Name: "test_MC_4p_noisy_single.m"
% Info: Script to run the Monte Carlo Simulation 
%       of plasma parameters estimation with the Levenberg-Marquardt 
%       optimization algorithm
% -------------------------------------------------------------------------
% "Levenberg-Marquardt optimization algorithm implementation to obtain the 
% Ionosphere parameters from the Incoherent Scatter Radar (ISR)" 
% -------------------------------------------------------------------------
% Based on:
%   Henri Gavin, Dept. Civil & Environ. Engineering, Duke Univ. 22 Sep 2013
%   http://www.duke.edu/~hpgavin/lm.m
%   http://people.duke.edu/~hpgavin/ce281/lm.pdf
%
% -------------------------------------------------------------------------
% Author: Miguel Martínez Ledesma (University of Chile)  
% Email: miguel.martinez@ing.uchile.cl
% Date: 19 November 2017
% -------------------------------------------------------------------------


%% INITIALIZATION

% Clean Workspace and Screen.-
clear all
clc

%% LIBRARIES

% Incoherent Scatter Radar Model Path.-
addpath('Model_ISR')
% LM Code Path.-
addpath('LM')

%% CONFIGURATIONS

% Objective function constants.-
function_options={};
function_options.radarFreq = 450e6;    % Radar frequency in Hz
function_options.NumPoints = 50;     % Length desired of the output vectors
function_options.maxFreq   = 10e3;     % Maximum frequency desired of the output vector

%% GENERATE DATA TO FIT

% Data generate.-
%
% True parameter values
%          [ n                Te    Ti   p ]
par_true = [ log10(2.7e+11)  520  480 0.51 ]; 
  
function_options.Vi=320;
 
% Function call
[input_Spectrum, frequency]= lm_func([], par_true', function_options);

%% Generate Spectrum With Noise.-

FluctuationPercentage=0.1;%10;%75%50%25%10;

standarddeviation=(FluctuationPercentage/100)*max(abs(input_Spectrum));

%Fixed fuctiation with frequency
noise = standarddeviation.*randn(length(input_Spectrum),1);

% Add Noise
noise_Spectrum=(input_Spectrum + noise);

%% INITIAL SEARCH PARAMETERS

initialSearch_parameters = [ log10(1e09)  300  300 0 ]; % initial guess for parameter values

%% L-M ALGORITHM OPTIONS

% Load LM Configuration.-
lm_configuration_parameters_4p

LM_configuration.options(1)=0; %Graphic output
LM_configuration.options(2)=100;%500; %Number of iterations (original is 2000)
LM_configuration.weight=standarddeviation^2; %Chi-squared Weigthing

%% SEPARATED LEVENBERG-MARQUARDT OPTIMIZATION ALGORITHM
[ fit_parameters, fit_value, sigma_parameters, sigma_Spectrum, correlation, R_2, convergence_history, secondselapsed, flag_stop, covar] = ...
            lm ('lm_func', initialSearch_parameters', frequency, noise_Spectrum, ...
            LM_configuration.weight, LM_configuration.fractional_increment_p, ...
            LM_configuration.par_min, LM_configuration.par_max, ...
            function_options, LM_configuration.options);


fit_Spectrum = lm_func(frequency,fit_parameters,function_options);
        
%% Change Electron density scale

par_true(1)=10.^par_true(1);
fit_parameters=fit_parameters';
fit_parameters(:,1)=10.^fit_parameters(:,1);

%% SHOW GRAPHIC

figure
hold on
plot(frequency./1e3,abs(input_Spectrum),'-b','linewidth',2)
plot(frequency./1e3,abs(noise_Spectrum),'-k')
plot(frequency./1e3,abs(fit_Spectrum),'-r','linewidth',2)
xlabel('Frequency [kHz]')
ylabel('Spectrum Magnitude')
box on
legend({'Original' 'Noisy' 'Fitted'})

%% PRINT RESULT

parameter_error=100*(fit_parameters'-par_true')./par_true';
parameter_string={'Ne'; 'Te'; 'Ti'; 'p '};
fprintf('\n');
fprintf('Fluctuation Percentage: %.1f %%\n',FluctuationPercentage)
for i=1:length(par_true)
    fprintf('Parameter %s: \tOriginal: %10.4f - \tFitting: %10.4f - \tError: %10.5f %%\n',parameter_string{i},par_true(i),fit_parameters(i),parameter_error(i))
end

%% CLOSE LIBRARIES

% Remove Paths.-
rmpath('Model_ISR')
rmpath('LM')
