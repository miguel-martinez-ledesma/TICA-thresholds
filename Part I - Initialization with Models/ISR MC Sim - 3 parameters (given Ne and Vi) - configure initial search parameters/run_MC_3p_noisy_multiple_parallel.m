
%% 
% Name: "run_MC_3p_noisy_multiple_parallel.m"
% Info: Script to run the Monte Carlo Simulation 
%       of plasma parameters estimation with the Levenberg-Marquardt 
%       optimization algorithm multiple times with different noise
%       executing it in parallel
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
%clc

%% PATHS
% Incoherent Scatter Radar Model Path.-
addpath('Model_ISR')
% LM Code Path.-
addpath('LM')

%% Parallelization.-

 Parallel_on  = 1; %Paralelization ON.- 
% Parallel_on = 0; %Paralelization OFF.-    

Num_Parallel_Executions = 5; % Number of Parallel Executions

if Parallel_on == 1
    %distcomp.feature( 'LocalUseMpiexec', false )
    
    desired_Cluster = parcluster('local');
    desired_Cluster.NumWorkers = Num_Parallel_Executions;  
    saveProfile(desired_Cluster);                          
    
    current_Cluster = matlabpool('size');
    
    if current_Cluster  ~= desired_Cluster.NumWorkers
        if current_Cluster > 0
            matlabpool('close');
        end
        %parpool('local',myCluster.NumWorkers)  %Matlab 2013b and further
        matlabpool('open',desired_Cluster.NumWorkers); %Matlab 2013a
    end
else
    if matlabpool('size')>0
        matlabpool('close');
    end
end


%% CONFIGURATIONS

% Number of Plasma parameters to fit
Number_PlasmaParameters=3;

% Noise Configuration
Number_Fluctuations=[100 50 20 10 5 2 1 0.5 0.2 0.1 0.05 0.02 0.01 0.001];

% Initial parameters range to true value
Number_Initial_Random_Range_Percentages=[100 50 40 30 20 10 5 1]/100;

% Define Number of parameter sets to evaluate.-
Number_InputParameters = 1000;

% Simulation Number of repetitions with randomly distributed start points  
Number_Random_Realizations = 250;

%% ISR RADAR CONFIGURATIONS

% Objective function options.-
function_options={};
function_options.radarFreq = 450e6;     % Radar frequency in Hz
function_options.NumPoints = 50;        % Length desired of the output vectors
function_options.maxFreq   = 10e3;      % Maximum frequency desired of the output vector

%% LM Configurations

% Load LM Configuration.-
lm_configuration_parameters_3p

LM_configuration.options(2)=100; %Maximum Number of Iterations 

LM_configuration.options(5)=100/(function_options.NumPoints+Number_PlasmaParameters+1); %Chi^2 minimum with normalization 

%% Loop for Initial Ranges

for ini_range=1:length(Number_Initial_Random_Range_Percentages)
    initial_random_range_percentage=Number_Initial_Random_Range_Percentages(ini_range)
    
    %% Loop For Different Fluctuation Numbers

    for fluctnum=1:length(Number_Fluctuations)
        FluctuationPercentage=Number_Fluctuations(fluctnum)


        %% Random Number Generator Seed Reset.-
        rng('default')

        %% Determine the Input Parameters to Evaluate.-

        % Load Input True Plasma Parameters.-
        load('True_Rand_Values.mat')

        %% Initialize Storage Struct

        LM_results={};

        %%    
        for j=1:Number_InputParameters

            % Initialize Timer.-
            timer_execution=tic();

            % Print Execution Number.-
            fprintf('\n \t Execution: %i\n',j)       

            %% Generate Struct with Input Plasma Parameters.-
            true_input_parameters(1)=True_Parameters.Te(j);
            true_input_parameters(2)=True_Parameters.Ti(j);
            true_input_parameters(3)=True_Parameters.p(j);

            function_options.Ne=10^11;%True_Parameters.Ne(j);
            function_options.Vi=0;%True_Parameters.Vi(j);

            %% Random Set of initial parameter.-

            initialSearch_parameters = [ ...
                true_input_parameters(1) + initial_random_range_percentage*((LM_configuration.par_max(1)-LM_configuration.par_min(1))*(rand(1,Number_Random_Realizations)-0.5));...
                true_input_parameters(2) + initial_random_range_percentage*((LM_configuration.par_max(2)-LM_configuration.par_min(2))*(rand(1,Number_Random_Realizations)-0.5));... 
                true_input_parameters(3) + initial_random_range_percentage*((LM_configuration.par_max(3)-LM_configuration.par_min(3))*(rand(1,Number_Random_Realizations)-0.5))]; 

            %% Generate Spectrum from input plasma parameters.-
            [input_Spectrum, frequency]= lm_func([], true_input_parameters', function_options);

            %% Generate Spectrum With Noise.-
            standarddeviation=(FluctuationPercentage/100)*max(abs(input_Spectrum));

            if standarddeviation==0
                LM_configuration.weight=1;
            else
                LM_configuration.weight=1/(standarddeviation); %Chi-squared Weigthing
            end

            %Fixed fuctiation with frequency
            noise = standarddeviation.*randn(length(input_Spectrum),Number_Random_Realizations);

            %Generate Noisy Spectrum
            noisy_Spectrum=[];
            for realization = 1:Number_Random_Realizations
                noisy_Spectrum(:,realization) = input_Spectrum + noise(:,realization);
            end

            %% Reserve Memory.-
            save_fit        = [];
            save_par        = [];
            save_iterations = [];
            save_flag_stop  = [];
            save_execution_time = [];

            saveall_fit = [];
            saveall_par = [];
            saveall_iterations = [];
            saveall_flag_stop = [];
            saveall_execution_time = [];

            %% MAIN OPTIMIZATION (Separated L-M).-

            parfor realization = 1:Number_Random_Realizations
            %for realization = 1:Number_Random_Realizations 

                [ fit_parameters, fit_value, sigma_parameters, sigma_Spectrum, correlation, R_2, convergence_history, secondselapsed, flag_stop, covar] = ...
                    lm ('lm_func', initialSearch_parameters(:,realization)', frequency, noisy_Spectrum(:,realization), ...
                    LM_configuration.weight, LM_configuration.fractional_increment_p, ...
                    LM_configuration.par_min, LM_configuration.par_max, ...
                    function_options, LM_configuration.options);

                % flag_stop  : Information of Stop Condition.-
                %             (1) Convergence in Gradient: epsilon1.-
                %             (2) Convergence in Parameters: epsilon2.-
                %             (3) Convergence in Chi-square: epsilon3.-
                %             (4) Maximum Number of Iterations Reached.- 
                %             (5) Floating Point Error (Infinite difference).- 
                %             (6) Initial Guess is Extremely Close to Optimal.- 
                %             (7) Maximum Number of Unaccepted Iterations.- 

                % Save Iteration Number.-
                save_fit(realization)               = fit_value;
                save_par(realization,:)             = fit_parameters;
                save_iterations(realization)        = length(convergence_history);
                save_flag_stop(realization)         = flag_stop;
                save_execution_time(realization)    = secondselapsed;
                save_sigma_parameters(realization,:)= sigma_parameters;
                save_R_2(realization)               = R_2;
                save_covar(realization,:,:)         = covar;

            end

            %% SAVE
            % Save Realization in Output Struct.-
            LM_results.fit_value(j,1:Number_Random_Realizations)         = save_fit(1:Number_Random_Realizations);
            LM_results.fit_parameters(j,1:Number_Random_Realizations,1:Number_PlasmaParameters) ...
                                                                         = save_par(1:Number_Random_Realizations,1:Number_PlasmaParameters);
            LM_results.iterations(j,1:Number_Random_Realizations)        = save_iterations(1:Number_Random_Realizations);
            LM_results.flag_stop(j,1:Number_Random_Realizations)         = save_flag_stop(1:Number_Random_Realizations);
            LM_results.execution_time(j,1:Number_Random_Realizations)    = save_execution_time(1:Number_Random_Realizations);

            LM_results.sigma_parameters(j,1:Number_Random_Realizations,1:Number_PlasmaParameters) ...
                                                                        = save_sigma_parameters(1:Number_Random_Realizations,1:Number_PlasmaParameters);
            LM_results.R_2(j,1:Number_Random_Realizations)              = save_R_2(1:Number_Random_Realizations);
            LM_results.covar(j,1:Number_Random_Realizations,:,:)        = save_covar(1:Number_Random_Realizations,:,:);

            % All Realizations Execution Time.-
            realization_time = toc(timer_execution);
            LM_results.realization_time(j) = realization_time;
            fprintf('Elapsed time is %.2f seconds \n',realization_time);

            % Print Optimum Values.-
            fprintf('LM Fit Value (%d) = [ ',j)
            for execution = 1:1:Number_Random_Realizations 
                fprintf('%.3d - ',LM_results.fit_value(j,execution))
            end
            fprintf(']')
            fprintf('\n')

            %% SAVE RESULT
            if mod(j,200)==1
                save(sprintf('%ip_%.1ernd_%irep_%iparam_%.2fini.mat',Number_PlasmaParameters,FluctuationPercentage,Number_Random_Realizations,Number_InputParameters,initial_random_range_percentage),...
                    'Number_InputParameters','Number_Random_Realizations','LM_results','LM_configuration','initialSearch_parameters');
            end
        end
        fprintf('DONE!\n')

        %% SAVE RESULT
        save(sprintf('%ip_%.1ernd_%irep_%iparam_%.2fini.mat',Number_PlasmaParameters,FluctuationPercentage,Number_Random_Realizations,Number_InputParameters,initial_random_range_percentage),...
            'Number_InputParameters','Number_Random_Realizations','LM_results','LM_configuration','initialSearch_parameters');

    end
    
end    
fprintf('DONE ALL!!!\n')

%% Finalize Execution and Close Everything

% Paralelization Tool (Close Matlabpool).-
if Parallel_on == 1
    % delete(gcp)
    matlabpool close
end

% Remove Paths.-
rmpath('Model_ISR')
rmpath('LM')