
%% 
% Name: "test_MC_4p_noisy_multiple.m"
% Info: Script to run the Monte Carlo Simulation 
%       of plasma parameters estimation with the Levenberg-Marquardt 
%       optimization algorithm multiple times with different noise
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

clear all
clc

%% LIBRARIES

% Incoherent Scatter Radar Model Path.-
addpath('Model_ISR')
% LM Code Path.-
addpath('LM')

%% CONFIGURATIONS

% Noise Configuration
FluctuationPercentage=10;

% Simulation Number of repetitions
Num_Iterations=250;

% Objective function constants.-
function_options={};
function_options.radarFreq = 450e6;    % Radar frequency in Hz
function_options.NumPoints = 50;       % Length desired of the output vectors
function_options.maxFreq   = 10e3;     % Maximum frequency desired of the output vector

%% GENERATE DATA TO FIT

% Data generate.-
%
% True parameter values
%          [ n                Te    Ti   p ]
par_true = [ log10(2.7e+11)  540  480 0.51 ]; 
  
function_options.Vi=320;
 
% Function call
[input_Spectrum, frequency]= lm_func([], par_true', function_options);

%% Generate Spectrum With Noise.-

standarddeviation=(FluctuationPercentage/100)*max(abs(input_Spectrum));

%Fixed fuctiation with frequency
noise = standarddeviation.*randn(Num_Iterations,length(input_Spectrum));

%% INITIAL SEARCH PARAMETERS

initialSearch_parameters = [ log10(1e09) 300  300 0 ]; % initial guess for parameter values

%% L-M ALGORITHM OPTIONS

% Load LM Configuration.-
lm_configuration_parameters_4p

LM_configuration.options(1)=0; %Graphic output
LM_configuration.options(2)=100; %Number of iterations 

% Weighting vector for least squares fit ( weight >= 0 ) ...
if standarddeviation==0
    LM_configuration.weight=1;
else
    LM_configuration.weight=1/(standarddeviation); %Standard Deviation Weigthing, later Chi-Squared
end        

%% SEPARATED LEVENBERG-MARQUARDT OPTIMIZATION ALGORITHM
noise_Spectrum=[];
fit_Spectrum=[];
fit_parameters=[];
sigma_parameters=[];
sigma_Spectrum=[];
correlation={};
flag_stop=[];
secondselapsed=[];
iterations=[];

for i=1:Num_Iterations
    % Add Noise
    noise_Spectrum(i,:)=(input_Spectrum + noise(i,:)');

    % Fit Parameters
    [ fit_parameters(i,:), fit_value, sigma_parameters(i,:), sigma_Spectrum(i,:), correlation{i}, R_2, convergence_history, secondselapsed(i), flag_stop(i), covar] = ...
        lm ('lm_func', initialSearch_parameters', frequency, noise_Spectrum(i,:), ...
        LM_configuration.weight, LM_configuration.fractional_increment_p, ...
        LM_configuration.par_min, LM_configuration.par_max, ...
        function_options, LM_configuration.options);

    % Generate Fitted Spectrum
    fit_Spectrum(i,:) = lm_func(frequency,fit_parameters(i,:)',function_options);
end

%% Change Electron density scale

par_true(1)=10.^par_true(1);
fit_parameters(:,1)=10.^fit_parameters(:,1);

%% SHOW SPECTRUM GRAPHIC 

figure
hold on
for i=1:Num_Iterations
    plot(frequency./1e3,abs(noise_Spectrum(i,:)),'-c','linewidth',1)
    plot(frequency./1e3,abs(fit_Spectrum(i,:)),'-g','linewidth',1)
end
h1=plot(frequency./1e3,abs(input_Spectrum),'-b','linewidth',2);
h2=plot(frequency./1e3,abs(mean(noise_Spectrum)),'-k','linewidth',2);
h3=plot(frequency./1e3,abs(mean(noise_Spectrum))+2*abs(std(noise_Spectrum)),':k','linewidth',2);
h4=plot(frequency./1e3,abs(mean(fit_Spectrum)),'-r','linewidth',2);
h5=plot(frequency./1e3,abs(mean(fit_Spectrum))+2*abs(std(fit_Spectrum)),':r','linewidth',2);
plot(frequency./1e3,abs(mean(noise_Spectrum))-2*abs(std(noise_Spectrum)),':k','linewidth',2)
plot(frequency./1e3,abs(mean(fit_Spectrum))-2*abs(std(fit_Spectrum)),':r','linewidth',2)
plot(frequency./1e3,abs(mean(noise_Spectrum)),'-k','linewidth',2)
xlabel('Frequency [kHz]')
ylabel('ISR Spectrum Magnitude')
ylim([0 max(max(noise_Spectrum))])
box on
legend([h1 h2 h3 h4 h5],{'Original Spectrum' 'Noisy Spectrum Mean (\mu_{noise})' 'Noisy Spectrum Std. Dev. (\sigma_{noise})' 'Fitted Spectrum Mean (\mu_{fitted})' 'Fitted Spectrum Std. Dev. (\sigma_{fitted})'})
title(sprintf('Noisy Spectrum Fitting (%i noisy samples and %.1e%% fluctuation)',Num_Iterations,FluctuationPercentage))

%% PRINT STATISTICS

for i=1:Num_Iterations
    parameter_error(i,:)=100*(fit_parameters(i,:)'-par_true')./par_true';
end
parameter_error_mean=mean(parameter_error);
parameter_error_std=std(parameter_error);

fit_parameters_mean=mean(fit_parameters);
fit_parameters_std=std(fit_parameters);

parameter_string={'Ne'; 'Te'; 'Ti'; 'p '};
fprintf('\n');
fprintf('Fluctuation Percentage: %.4f %%\n',FluctuationPercentage)
fprintf('Number of Noisy Samples: %i\n',Num_Iterations)
for i=1:length(par_true)
    fprintf('Parameter %s: \tOriginal: %10.4f \t- Fitting: Mean=%10.4f Std=%8.4f \t- Error: Mean=%8.5f%% Std=%8.5f\n',...
        parameter_string{i},par_true(i),fit_parameters_mean(i),fit_parameters_std(i),parameter_error_mean(i),parameter_error_std(i))
end
fprintf('Fit Value: \tMean=%.6e \tStd=%.6e\n',mean(fit_value),std(fit_value));

%% SHOW HISTOGRAM

parameter_string={'Ne'; 'Te'; 'Ti'; 'p '};
visual_added=[25 25 25 0.1 0.1];
hist_num=50;
figure

for i=1:length(par_true)

    subplot(length(par_true),1,i)
    hold on
    parameter=i;
    [f,binCtrs]=hist(fit_parameters(:,parameter),hist_num);
    bar(binCtrs,f/max(f)); 
    maxhist=max(f/trapz(binCtrs,f));
    plot([par_true(parameter) par_true(parameter)],[0 1],'-r','linewidth',2)
    plot([fit_parameters_mean(parameter) fit_parameters_mean(parameter)],[0 1],'--b','linewidth',2)
    gauss_space=linspace(min(fit_parameters(:,parameter))-visual_added(i),max(fit_parameters(:,parameter))+visual_added(i),100);
    gauss_value=normpdf(gauss_space,fit_parameters_mean(parameter),fit_parameters_std(parameter));
    plot(gauss_space,gauss_value.*(1/max(gauss_value)),'--b','linewidth',1)
    box on
    axis tight
    xlabel(parameter_string{i})
    ylabel('Normalized Histogram')
    legend({sprintf('Histogram (%i samples)',Num_Iterations) ...
            sprintf('True Value: %10.4e', par_true(parameter))...
            sprintf('Fitted: \\mu=%.5e \\sigma=%.5e',fit_parameters_mean(parameter),fit_parameters_std(parameter))})
    ylim([0 1])
end
subplot(length(par_true),1,1)
title(sprintf('Results Fitting with Levenberg-Marquardt - Noise fluctuation %.4f%%',FluctuationPercentage))


%% CLOSE LIBRARIES

% Remove Paths.-
rmpath('Model_ISR')
rmpath('LM')
