% statistical_analysis_function.m
%
% Function to analyze Monte Carlo results of an ISR estimation simulation
% -------------------------------------------------------------------------
% Inputs: 
%   - filename: MonteCarlo simulation file
%   - FluctuationPercentage: Signal Noise of simulation
%   - PlasmaParameters: Name of Plasma Parameters simulated 
%       example: {'Ne','Te','Ti','Vi','p'}
%   - Selection_Parameter: Name of parameter to determine correctness 
%       example: 'p'
%
% -------------------------------------------------------------------------
% Author: Miguel Martínez Ledesma (University of Chile)  
% Email: miguel.martinez@ing.uchile.cl
% Date: 24 November 2016
% -------------------------------------------------------------------------

function [  Number_InputParameters,...
            LM_results, ...
            LM_configuration, ...
            LM_results_valid, ...
            LM_results_correct] = statistical_analysis_function(filename,FluctuationPercentage,PlasmaParameters, Selection_Parameter)

%% Define Chi-square limit (for valid convergence selection)

X2_min=100;

%% Verify if input file exists            
if exist(filename, 'file') == 2
    load(filename);
else
    error ('Error: file "%s" does not exist!!\n',filename)
end

%% Add Info to the Results Stuct 

load('True_Rand_Values.mat')

Number_PlasmaParameters=length(PlasmaParameters);

%Load Number of Stored Parameters
[Number_InputParameters, Number_Realizations]=size(LM_results.fit_value);

if Number_Realizations~=Number_Random_Realizations;
    Number_InputParameters=Number_Realizations; 
end

if Number_InputParameters<=10
    warning ('Warning: file "%s" contain less than 10 values\n',filename)
    LM_configuration=[];
    LM_results_valid=[];
    LM_results_correct=[];
    return
end

%% Input Parameters Selection

ListOfPlasmaParameters={'Ne','Te','Ti','Vi','p'};
ListOfPlasmaParametersIndex=[NaN,NaN,NaN,NaN,NaN];
Selection_Parameter_index=NaN;

for i=1:length(PlasmaParameters)
    for j=1:length(ListOfPlasmaParameters)
        if strcmp(PlasmaParameters{i},ListOfPlasmaParameters{j})==1
            ListOfPlasmaParametersIndex(j)=i;
        end
        if strcmp(PlasmaParameters{i},Selection_Parameter)==1
            Selection_Parameter_index=i;
        end
    end
end

if isnan(Selection_Parameter_index)
    warning ('Warning: Wrong Name for Selection Parameter\n')
    LM_configuration=[];
    LM_results_valid=[];
    LM_results_correct=[];
    return
end

True_Parameters_Vector=[];
if ~isnan(ListOfPlasmaParametersIndex(1))
    True_Parameters_Vector=[True_Parameters_Vector ...
                        True_Parameters.Ne(1:Number_InputParameters) ];
    LM_results.fit_parameters(:,:,1)=10.^LM_results.fit_parameters(:,:,1);
end
if ~isnan(ListOfPlasmaParametersIndex(2))
    True_Parameters_Vector=[True_Parameters_Vector ...
                        True_Parameters.Te(1:Number_InputParameters) ];
end
if ~isnan(ListOfPlasmaParametersIndex(3))
    True_Parameters_Vector=[True_Parameters_Vector ...
                        True_Parameters.Ti(1:Number_InputParameters) ];
end
if ~isnan(ListOfPlasmaParametersIndex(4))
    True_Parameters_Vector=[True_Parameters_Vector ...
                        True_Parameters.Vi(1:Number_InputParameters) ];
end
if ~isnan(ListOfPlasmaParametersIndex(5))
    True_Parameters_Vector=[True_Parameters_Vector ...
                        True_Parameters.p(1:Number_InputParameters) ];
end

LM_results.true_parameters=reshape(repmat(True_Parameters_Vector,Number_Random_Realizations,1),[Number_InputParameters Number_Random_Realizations Number_PlasmaParameters]);
LM_results.err_parameters=(LM_results.true_parameters)-(LM_results.fit_parameters);
%Delete covar matrix from results
LM_results=rmfield(LM_results,'covar');

clear True_Parameters_Vector True_Parameters ListOfPlasmaParameters ListOfPlasmaParametersIndex
clear filename


%% Filter Chi-square Minimum Limit: VALID results

LM_results_valid=LM_results;
for i=1:Number_InputParameters
    min_index=find(LM_results_valid.fit_value(i,:)>=X2_min);
    LM_results_valid.fit_value(i,min_index)=NaN;
    LM_results_valid.fit_parameters(i,min_index,1:Number_PlasmaParameters)=NaN;
    LM_results_valid.sigma_parameters(i,min_index,1:Number_PlasmaParameters)=NaN;
    LM_results_valid.true_parameters(i,min_index,1:Number_PlasmaParameters)=NaN;
    LM_results_valid.err_parameters(i,min_index,1:Number_PlasmaParameters)=NaN;
    LM_results_valid.flag_stop(i,min_index)=NaN;
    LM_results_valid.execution_time(i,min_index)=NaN;
    LM_results_valid.R_2(i,min_index)=NaN;
end
clear min_index i

 
%% Calculate Mixture of 2 Gaussians using Expectation Maximization EM algorithm

X=LM_results.true_parameters(:,:,Selection_Parameter_index);
Y=LM_results_valid.err_parameters(:,:,Selection_Parameter_index);

%Sort by Order
[~,I]=sort(X(:,1));

X=X(I,:);
Y=Y(I,:);

Ysize=1000;

Y__=linspace(-1,1,Ysize);

for j=1:length(X(:,1))

   %delete peaks at limit values(fitting at 0 or 1) at original data
   Y_filtered=Y(j,:)';
    if FluctuationPercentage>0.05 
        min_X=X(j,1)-0.005;
        max_X=X(j,1)-1+0.005;
        Y_filtered(Y_filtered>=min_X)=NaN;
        Y_filtered(Y_filtered<=max_X)=NaN;
    end

    %calculate histogram
    [f]=hist(Y_filtered,Y__);
    f=(f./sum(f));

    % Use Expectation-Maximization algorithm to compute PDF
    iterations=500;
    num_gaussians=2;
    log_likelihood_step=1e-30;

    initial_guess={};

    FluctuationList=[100 50 20 10 5 2 1 0.5 0.2 0.1 0.05 0.02 0.01 0.001];
    k=find(FluctuationList==FluctuationPercentage);
    if ~isempty(k)
        p_est_guess = 0.5+0.45*k/length(FluctuationList);
    else
        p_est_guess = 0.5;
    end

    if length(Y_filtered(~isnan(Y_filtered)))>1
        
        initial_guess.W=[p_est_guess 1-p_est_guess];
        initial_guess.M=[0 2*X(j,1)-1];
        initial_guess.V=zeros(1,1,2);
        initial_guess.V(:,:,1:2)=FluctuationPercentage/10;

        [weight,mu,variance,L] = EM_GM_fast(Y_filtered(~isnan(Y_filtered)),num_gaussians,log_likelihood_step,iterations,0,initial_guess);
        sigma=sqrt(reshape(variance,[1 2]));

        %find points with same mean and find only one gaussian
        if abs(mu(1)-mu(2))<0.1 || isnan(mu(2)) || weight(1)>0.975

            num_gaussians=1;
            [weight,mu,variance,L] = EM_GM_fast(Y_filtered(~isnan(Y_filtered)),num_gaussians,log_likelihood_step,iterations,0,initial_guess);
            sigma=sqrt(variance);

            mu = [mu(1) NaN];
            weight = [1 0];
            sigma = [sigma(1) NaN];

        end
    else
        mu = [NaN NaN];
        weight = [1 0];
        sigma = [NaN NaN];
    end

    %determine if gaussians is 0 or 1 group
    if isnan(mu(2)) && abs(mu(1))> 0.5 %if only one gaussian with mean far from 0
            mu=[NaN mu(1)];
            weight=[NaN weight(1)];
            sigma=[NaN sigma(1)];
    elseif abs(mu(1))>abs(mu(2)) %select the gaussian with smaller absolute mean
        mu=mu([2 1]);
        weight=weight([2 1]);
        sigma=sigma([2 1]);
    end



    estimated_alpha(j)=weight(1);
    estimated_mean_0(j)=mu(1);
    estimated_mean_1(j)=mu(2);
    estimated_sigma_0(j)=sigma(1);
    estimated_sigma_1(j)=sigma(2);
end
        

%% Cluster Correct and Incorrect Data: CORRECT results
    
LM_results_correct=LM_results_valid;

% X: true parameters, Y: estimated parameters
X=LM_results.true_parameters(:,:,Selection_Parameter_index);
Y=LM_results_valid.fit_parameters(:,:,Selection_Parameter_index);

%Sort X by Order
[~,I]=sort(X(:,1));
X=X(I,:);
Y=Y(I,:);

for j=1:length(X(:,1))

    if ~isnan(estimated_sigma_1(j)) && ~isnan(estimated_mean_1(j))
        if estimated_alpha(j)<1 && (estimated_sigma_1(j)~=0)
            sigma=zeros(1,1,2);
            sigma(:,:,1)=estimated_sigma_0(j)^2;
            sigma(:,:,2)=estimated_sigma_1(j)^2;
            gm_obj = gmdistribution([X(j)-estimated_mean_0(j) X(j)-estimated_mean_1(j)]',sigma,[estimated_alpha(j) 1-estimated_alpha(j)]);
            idx = cluster(gm_obj,Y(j,:)');

            correct_index=find(idx==1);
            incorrect_index=find(idx==2);

           
            LM_results_correct.fit_value(I(j),incorrect_index)=NaN;
            LM_results_correct.fit_parameters(I(j),incorrect_index,1:Number_PlasmaParameters)=NaN;
            LM_results_correct.sigma_parameters(I(j),incorrect_index,1:Number_PlasmaParameters)=NaN;
            LM_results_correct.true_parameters(I(j),incorrect_index,1:Number_PlasmaParameters)=NaN;
            LM_results_correct.err_parameters(I(j),incorrect_index,1:Number_PlasmaParameters)=NaN;
            LM_results_correct.flag_stop(I(j),incorrect_index)=NaN;
            LM_results_correct.execution_time(I(j),incorrect_index)=NaN;
            LM_results_correct.R_2(I(j),incorrect_index)=NaN;

        end
    end
    
end
    
