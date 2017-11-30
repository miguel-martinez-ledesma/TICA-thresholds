% create_probability_curves_3p.m
%
% Script to analyze Monte Carlo results of the ISR estimation simulation
% to compare the effect of wrong initializations of optimization parameters 
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

Number_InitialPercentage_List=[100 50 40 30 20 10 5 1]./100;
FluctuationList=[100 10 1 0.1 0.01 0.001];

for n=1:length(Number_InitialPercentage_List)
    initial_random_range_percentage=Number_InitialPercentage_List(n);

    for f=1:length(FluctuationList)
        FluctuationPercentage=FluctuationList(f);

		
        PlasmaParameters={'Te','Ti','p'};
        Selection_Parameter='p';

        Number_Random_Realizations = 250;
        Number_InputParameters = 1000;
        Number_PlasmaParameters = length(PlasmaParameters);
        
        filename=sprintf('%ip_%.1ernd_%irep_%iparam_%.2fini.mat',Number_PlasmaParameters,FluctuationPercentage,Number_Random_Realizations,Number_InputParameters,initial_random_range_percentage);
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
                warning ('file "%s" does not contain full dataset (50%% less)\n',filename)
                wrong_file_content=1;
            end
       	else
            warning ('file "%s" does not exist\n',filename)
            wrong_file_content=1;
        end
		
		if wrong_file_content==1
			%When an error has been found reading the file
			
            Solved_all(n,f)=NaN;
            Correct_all(n,f)=NaN;
            Iterations_all(n,f)=NaN;
            
            Solved(n,f,1:Number_InputParameters)=nan(Number_InputParameters,1);
            Correct(n,f,1:Number_InputParameters)=nan(Number_InputParameters,1);
            Iterations(n,f,1:Number_InputParameters)=nan(Number_InputParameters,1);
        
        else
			%Save results in an array for comparison
			
			Solved(n,f,1:Number_InputParameters)=sum(~isnan(LM_results_valid.fit_value),2);
			Correct(n,f,1:Number_InputParameters)=sum(~isnan(LM_results_correct.fit_value),2);
			Iterations(n,f,1:Number_InputParameters)=sum(LM_results.iterations,2);
			
			Solved_all(n,f)=length(find(~isnan(LM_results_valid.fit_value)));
			Correct_all(n,f)=length(find(~isnan(LM_results_correct.fit_value)));
			Iterations_all(n,f)=sum(sum(LM_results.iterations,2));
			
		end
		clear filename

    end
end
clear f n wrong_file_content


%% Save all results

save('Analyzed_Results.mat',...
    'Solved_all','Correct_all','Iterations_all',...
    'Solved','Correct','Iterations',...
    'PlasmaParameters','Number_PlasmaParameters','Selection_Parameter',...
    'FluctuationList','Number_Random_Realizations','Number_InputParameters',...
    'Number_InitialPercentage_List')

%% Remove Paths.-
rmpath('EM')

%% Load all results from file and clear everything else

clear all
load('Analyzed_Results.mat')

%Create Legend
legendstring={};
for n=1:length(Number_InitialPercentage_List)
    legendstring{n}=sprintf('Initial Range \\beta=%.0f %%',100*Number_InitialPercentage_List(n));
end
clear n

%Create colors for each Initial Search Range
newColors=jet(length(Number_InitialPercentage_List));


%% Show fluctuation results

figure
subplot(3,1,1)
set(gca,'colororder',newColors, 'NextPlot', 'replacechildren');
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
set(gca,'colororder',newColors, 'NextPlot', 'replacechildren');
hold on
plot(FluctuationList,100*Correct_all./(Solved_all),'o-','markersize',5,'linewidth',2)
grid on
box on
set(gca,'XScale','Log');
xlabel('Fluctuation Percentage ( \delta(%) )')
ylabel('P_{Correct} (%)')
ylim([40 100])

subplot(3,1,3)
set(gca,'colororder',newColors, 'NextPlot', 'replacechildren');
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
set(gca,'colororder',newColors, 'NextPlot', 'replacechildren');
hold on
plot(FluctuationList,Iterations_all/(Number_InputParameters*Number_Random_Realizations),'o-','linewidth',2)
grid on
box on
set(gca,'XScale','Log');
xlabel('Fluctuation Percentage ( \delta(%) )')
ylabel('Average Iterations (num)')
legend(legendstring,'Location','Best')

set(gcf,'Position', [20, 20, 600, 500])

