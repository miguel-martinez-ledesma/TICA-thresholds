% create_probability_curves_4p.m
%
% Script to analyze Monte Carlo results of the ISR estimation simulation
% given information of parameter Ne
% -------------------------------------------------------------------------
% Author: Miguel Martínez Ledesma (University of Chile)  
% Email: miguel.martinez@ing.uchile.cl
% Date: 24 November 2016
% -------------------------------------------------------------------------

%% Initialize
clear all

%% PATHS
% Expected Maximization Path.-
addpath('EM')

%% Initialize 
Solved=[];
Correct=[];
Iterations=[];
        
Solved_all=[];
Correct_all=[];
Iterations_all=[];

%% Load Results

FluctuationList=[100 50 20 10 5 2 1 0.5 0.2 0.1 0.05 0.02 0.01 0.001];

for f=1:length(FluctuationList)
    FluctuationPercentage=FluctuationList(f);


    PlasmaParameters={'Te','Ti','Vi','p'};
    Selection_Parameter='p';

    Number_Random_Realizations = 250;
    Number_InputParameters = 1000;
    Number_PlasmaParameters = length(PlasmaParameters);

    filename=sprintf('%ip_%.1ernd_%irep_%iparam.mat',Number_PlasmaParameters,FluctuationPercentage,Number_Random_Realizations,Number_InputParameters);
	filename

    %Load file with Monte Carlo results 
    wrong_file_content=0;
    if exist(filename, 'file') == 2

        [Number_InputParameters,...
        LM_results, ...
        LM_configuration, ...
        LM_results_valid, ...
        LM_results_correct] = statistical_analysis_function(filename,FluctuationPercentage,PlasmaParameters,Selection_Parameter);

        if length(LM_results.fit_value)<Number_InputParameters/2
            warning ('file "%s" does not contain full dataset (less than 50%%)\n',filename)
            wrong_file_content=1;
        end
    else
        warning ('file "%s" does not exist\n',filename)
        wrong_file_content=1;
    end

    if wrong_file_content==1
        %When an error has been found reading the file

        Solved_all(f)=NaN;
        Correct_all(f)=NaN;
        Iterations_all(f)=NaN;

        Solved(f,1:Number_InputParameters)=nan(Number_InputParameters,1);
        Correct(f,1:Number_InputParameters)=nan(Number_InputParameters,1);
        Iterations(f,1:Number_InputParameters)=nan(Number_InputParameters,1);

    else
        %Save results in an array for comparison

        Solved(f,1:Number_InputParameters)=sum(~isnan(LM_results_valid.fit_value),2);
        Correct(f,1:Number_InputParameters)=sum(~isnan(LM_results_correct.fit_value),2);
        Iterations(f,1:Number_InputParameters)=sum(LM_results.iterations,2);

        Solved_all(f)=length(find(~isnan(LM_results_valid.fit_value)));
        Correct_all(f)=length(find(~isnan(LM_results_correct.fit_value)));
        Iterations_all(f)=sum(sum(LM_results.iterations,2));

    end
    clear filename

end

clear f wrong_file_content


%% Save all results

save('Analyzed_Results_4p_given_Ne.mat',...
    'Solved_all','Correct_all','Iterations_all',...
    'Solved','Correct','Iterations',...
    'PlasmaParameters','Number_PlasmaParameters','Selection_Parameter',...
    'FluctuationList','Number_Random_Realizations','Number_InputParameters')

%% Remove Paths.-
rmpath('EM')

%% Load all results from file and clear everything else

clear all
load('Analyzed_Results_4p_given_Ne.mat')

%Create Legend
legendstring={};
legendstring{1}=sprintf('Analyzed Results 4p given N_e');


%% Show fluctuation results

figure
subplot(3,1,1)
hold on
plot(FluctuationList,100*Solved_all/(Number_InputParameters*Number_Random_Realizations),'o-','markersize',5,'linewidth',2)
grid on
box on
set(gca,'XScale','Log');
xlabel('Fluctuation Percentage ( \delta(%) )')
ylabel('P_{Valid} (%)')
legend(legendstring,'Location','Best')
ylim([40 100])

subplot(3,1,2)
hold on
plot(FluctuationList,100*Correct_all./(Solved_all),'o-','markersize',5,'linewidth',2)
grid on
box on
set(gca,'XScale','Log');
xlabel('Fluctuation Percentage ( \delta(%) )')
ylabel('P_{Correct} (%)')
ylim([40 100])

subplot(3,1,3)
hold on
plot(FluctuationList,100*Correct_all./(Number_InputParameters*Number_Random_Realizations),'o-','markersize',5,'linewidth',2)
grid on
box on
set(gca,'XScale','Log');
xlabel('Fluctuation Percentage ( \delta(%) )')
ylabel('P_{Valid & Correct} (%)')
ylim([40 100])

set(gcf,'Position', [20, 20, 650, 950])

%% Show Number of Iterations

figure
hold on
plot(FluctuationList,Iterations_all/(Number_InputParameters*Number_Random_Realizations),'o-','linewidth',2)
grid on
box on
set(gca,'XScale','Log');
xlabel('Fluctuation Percentage ( \delta(%) )')
ylabel('Average Iterations (num)')
legend(legendstring,'Location','Best')

set(gcf,'Position', [20, 20, 600, 500])

