% graphical_statistical_analysis_3p.m
%
% Function to analyze Monte Carlo results graphically
%
% -------------------------------------------------------------------------
% Author: Miguel Martínez Ledesma (University of Chile)  
% Email: miguel.martinez@ing.uchile.cl
% Date: 24 November 2016
% -------------------------------------------------------------------------

%% Initialize
clear all

%% Configuration
save_filename='Graphically_Analyzed_Results_3p_given_NeVi.mat';

PlasmaParameters={'Te','Ti','p'};
Selection_Parameter='p';

%% Load Results

Number_PlasmaParameters = length(PlasmaParameters);
Number_Random_Realizations = 250;
Number_InputParameters_Expected = 1000;
    
FluctuationList=[100 50 20 10 5 2 1 0.5 0.2 0.1 0.05 0.02 0.01 0.001];

Fluctuation_Colors=jet(length(FluctuationList));

%% Configurations

% Define Chi-square limit (for correct result selection)
X2_min=1e2;

%% Initialize Structures
Y_=[];
X_=[];
S1_mean=[];
S1_std=[];
S2_mean=[];
S2_std=[];

paramEsts=[];

estimated_alpha=[];
estimated_mean_1=[];
estimated_mean_0=[];
estimated_sigma_1=[];
estimated_sigma_0=[];
        
%% Loop for each result file
for k=1:length(FluctuationList)
    FluctuationPercentage=FluctuationList(k);

    %% Load File
    filename=sprintf('%ip_%.1ernd_%irep_%iparam_1segm0.10over.mat',Number_PlasmaParameters,FluctuationPercentage,Number_Random_Realizations,Number_InputParameters_Expected);
    filename
    if exist(filename, 'file') == 2
        load(filename);
    else
        warning ('Warning: file "%s" does not exist!!\n',filename)
        continue
    end
    clear filename

    %Reload Number of Stored Parameters
    Number_InputParameters=length(LM_results.fit_value);

    
    %% True Parameters Selection
    load('True_Rand_Values.mat')

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

    
    %% Filter Chi-square Minimum Limit

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
    
    %% Calculate Mixture of 2 Gaussians using Expectation Maximization algorithm
    
    %X: true parameters, Y: error of parameter estimations
    X=LM_results.true_parameters(:,:,Selection_Parameter_index);
    Y=LM_results_valid.err_parameters(:,:,Selection_Parameter_index);
    
    %Sort by Order
    [~,I]=sort(X(:,1));
    X=X(I,:);
    Y=Y(I,:);
    
    %histogram axis creation
    Ysize=1000;
    Y__=linspace(-1,1,Ysize);

    for j=1:length(X(:,1))
        
       %delete peaks at limit values(fitting at 0 or 1) from original data
       Y_filtered=Y(j,:)';
        if FluctuationPercentage>0.05 
            min_X=X(j,1)-0.005;
            max_X=X(j,1)-1+0.005;
            Y_filtered(Y_filtered>=min_X)=NaN;
            Y_filtered(Y_filtered<=max_X)=NaN;
        end

        %calculate histogram for graphical visualization
        [f]=hist(Y_filtered,Y__);
        f=(f./sum(f)); %Normalize histogram
        
        %% Use Expectation-Maximization algorithm to compute PDF
        if length(Y_filtered(~isnan(Y_filtered)))<=1
            mu=[NaN NaN];
            weight=[1 0];
            sigma=[NaN NaN];
        else
			%Configuration
			iterations=500;
			num_gaussians=2;
			log_likelihood_step=1e-30;

			%Initialization
			initial_guess={};
			p_est_guess=0.5+0.45*k/length(FluctuationList);
			initial_guess.W=[p_est_guess 1-p_est_guess];
			initial_guess.M=[0 2*X(j,1)-1];
			initial_guess.V=zeros(1,1,2);
			initial_guess.V(:,:,1:2)=FluctuationPercentage/10;
			
			%Expected Maximization
			[weight,mu,variance,L] = EM_GM_fast(Y_filtered(~isnan(Y_filtered)),num_gaussians,log_likelihood_step,iterations,0,initial_guess);
			sigma=sqrt(reshape(variance,[1 2]));
			
			%If gaussians have equal mean, then use only one gaussian
			if abs(mu(1)-mu(2))<0.1 || isnan(mu(2)) || weight(1)>0.975
				num_gaussians=1;
				[~,mu,variance,L] = EM_GM_fast(Y_filtered(~isnan(Y_filtered)),num_gaussians,log_likelihood_step,iterations,0,initial_guess);
				sigma=sqrt(variance);
			
				weight = [1 0];
				mu = [mu(1) NaN];
				sigma = [sigma(1) NaN];
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
		end
			
        %Save statistical params
        estimated_alpha(k,j)=weight(1);
        estimated_mean_0(k,j)=mu(1);
        estimated_mean_1(k,j)=mu(2);
        estimated_sigma_0(k,j)=sigma(1);
        estimated_sigma_1(k,j)=sigma(2);
        
        %% Plot Histogram and Estimated Statistical Parameters
        figure(1)
        
        %create PDF of correct and incorrect distributions together
        pdf_estimated = weight(1)*normpdf(Y__,mu(1),sigma(1));
        if ~isnan(weight(2)) && (weight(2)>0)
            pdf_estimated = pdf_estimated + weight(2)*normpdf(Y__,mu(2),sigma(2));
        end
        %Normalize PDF
        pdf_estimated = pdf_estimated/sum(pdf_estimated);
        
        %Clear screen
        hold off
        %Plot histogram
        plot(Y__,f,'-k')
        hold on
        
        %Show in red two gaussians and in green only one
        if ~isnan(weight(2)) && (weight(2)>0)
            plot(Y__,pdf_estimated,'-r','linewidth',2)
        else
            plot(Y__,pdf_estimated,'-g','linewidth',2)
        end
        
        %Show bars with mean values: blue correct, red incorrect
        plot([mu(1) mu(1)],[min(f) max(f)]./2,'--b','linewidth',2)
        plot([mu(2) mu(2)],[min(f) max(f)]./2,'--r','linewidth',2)
        
        %Show index of data analyzed (from -1 to 1)
        plot([2*X(j,1)-1 2*X(j,1)-1],[min(f) max(f)]./2,'-c','linewidth',2)
        
        %Show texts
        xlabel('Error in Ion Composition')
        ylabel('Histogram')
        if FluctuationPercentage>=1
            title(sprintf('Fluctuation Percentage \\delta=%i%%',FluctuationPercentage))
        elseif FluctuationPercentage>=0.01
            title(sprintf('Fluctuation Percentage \\delta=%.2f%%',FluctuationPercentage))
        else
            title(sprintf('Fluctuation Percentage \\delta=%.3f%%',FluctuationPercentage))
        end
        
        %Graphic limits and colors
        box on
        if max(f)>0
            ylim([0 max(f)])
        end
        set(gca,'Color',[1 1 1]);
        set(gcf,'Color','white')
        xlim([-1 1])
        
    end

end

%%    
save(save_filename,...
                    'PlasmaParameters','Selection_Parameter',...
                    'Number_InputParameters_Expected','Number_Random_Realizations','X2_min',...
                    'FluctuationList','Fluctuation_Colors',...
                    'estimated_alpha','estimated_mean_0','estimated_mean_1','estimated_sigma_0','estimated_sigma_1');
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
%% CLEAR UNUSED VARIABLES
%clear all
clearvars -except save_filename
load (sprintf(save_filename))                             

%% CLUSTER DATA


CountCorrect=zeros([1,length(FluctuationList)]);
CountIncorrect=zeros([1,length(FluctuationList)]);

for k=1:length(FluctuationList)
    FluctuationPercentage=FluctuationList(k);

    %% Load File
    Number_PlasmaParameters = length(PlasmaParameters);
    filename=sprintf('%ip_%.1ernd_%irep_%iparam_1segm0.10over.mat',Number_PlasmaParameters,FluctuationPercentage,Number_Random_Realizations,Number_InputParameters_Expected);
    filename
    if exist(filename, 'file') == 2
        load(filename);
    else
        warning ('Warning: file "%s" does not exist!!\n',filename)
        continue
    end
    clear filename

    %Reload Number of Stored Parameters
    Number_InputParameters=length(LM_results.fit_value);

    
    %% True Parameters Selection
    load('True_Rand_Values.mat')

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

    
    %% Filter Chi-square Minimum Limit

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
    
    %% Show Valid and Invalid Values
    
    X=LM_results.true_parameters(:,:,Selection_Parameter_index);
    Y=LM_results_valid.fit_parameters(:,:,Selection_Parameter_index);
    
    %Sort by Order
    [~,I]=sort(X(:,1));
    X=X(I,:);
    Y=Y(I,:);

    Y_correct=Y;
    Y_incorrect=Y;
    
    prob=estimated_alpha(k,:);
    M0=estimated_mean_0(k,:);
    M1=estimated_mean_1(k,:);
    S0=estimated_sigma_0(k,:);
	S1=estimated_sigma_1(k,:);
    
    Ylim=NaN(length(X(:,1)),2);
    
    for j=1:length(X(:,1))
        
        
        if ~isnan(S1(j)) && ~isnan(M1(j))
            if prob(j)<1 && (S1(j)~=0)
                sigma=zeros(1,1,2);
                sigma(:,:,1)=S0(j)^2;
                sigma(:,:,2)=S1(j)^2;
                gm_obj = gmdistribution([X(j)-M0(j) X(j)-M1(j)]',sigma,[prob(j) 1-prob(j)]);
                idx = cluster(gm_obj,Y(j,:)');

                correct_index=find(idx==1);
                incorrect_index=find(idx==2);
                
                Y_correct(j,incorrect_index)=NaN;
                Y_incorrect(j,correct_index)=NaN;

                CountCorrect(k)=CountCorrect(k)+length(find(~isnan(Y_correct(j,:))));
                CountIncorrect(k)=CountIncorrect(k)+length(find(~isnan(Y_incorrect(j,:))));
            else
                Y_incorrect(j,:)=NaN;
                CountCorrect(k)=CountCorrect(k)+length(find(~isnan(Y_correct(j,:))));
            end
        else
            Y_incorrect(j,:)=NaN;
            CountCorrect(k)=CountCorrect(k)+length(find(~isnan(Y_correct(j,:))));
        end
        
    end
    
    
    figure(2)
    subplot(ceil(length(FluctuationList)/2),2,k)
    hold on    
    
    plot(X,reshape(Y_incorrect,size(X)),'ob','markersize',1)
    plot(X,reshape(Y_correct,size(X)),'or','markersize',1)
    
    if FluctuationPercentage>=1
        title(sprintf('\\delta=%i%% (P_{correct}=%.3f%%)',FluctuationPercentage,100*CountCorrect(k)/(CountCorrect(k)+CountIncorrect(k))))
    elseif FluctuationPercentage>=0.01
        title(sprintf('\\delta=%.2f%%(P_{correct}=%.3f%%)',FluctuationPercentage,100*CountCorrect(k)/(CountCorrect(k)+CountIncorrect(k))))
    else
        title(sprintf('\\delta=%.3f%% (P_{correct}=%.3f%%)',FluctuationPercentage,100*CountCorrect(k)/(CountCorrect(k)+CountIncorrect(k))))
    end    
    
    box on
    grid on
    xlabel('True ion composition (p)')
    ylabel('Fitted p')
    ylim([0 1])
    set(gca,'Color',[1 1 1]);
    set(gcf,'Color','white')
    set(gcf,'Position', [100, 0, 800, 1300])
     
end

%%
figure(3)
hold on
plot(FluctuationList,100.*CountCorrect./(CountCorrect+CountIncorrect),'o-r','linewidth',2)
box on
%grid on
xlabel('Fluctuation Level Percentage (\delta)')
ylabel('Valid Selection of Ion Composition ( Percentage )')
set(gca,'XScale','Log');
set(gca,'Color',[1 1 1]);
set(gcf,'Color','white')

clear cbh liststr k prob S1 M1 S0 M0  X__ X Y i j Xsize divider valid_index invalid_index FluctuationPercentage
clear gm_obj idx Y_incorrect Y_correct sigma meanX index
clear CountIncorrect CountCorrect





%% SHOW GAUSSIAN STANDARD DEVIATION 

legendindex=[];
legendstring={};

figure(4)
hold on

S0_mean=nanmean(estimated_sigma_0,2)';
S1_mean=nanmean(estimated_sigma_1,2)';
for i=1:length(S0_mean)
    if (length(find(~isnan(estimated_sigma_0(i,:))))<=10)
       S0_mean(i)=NaN; 
    end
    if (length(find(~isnan(estimated_sigma_1(i,:))))<=10)
       S1_mean(i)=NaN; 
    end
end

%------------------------
h=plot(FluctuationList, S0_mean,'s-r','markersize',6,'linewidth',2);
legendindex(end+1)=h;
legendstring{end+1}='\sigma_{0}';

%------------------------
h=plot(FluctuationList(1:length(S1_mean)), S1_mean,'o-b','markersize',6,'linewidth',2);
legendindex(end+1)=h;
legendstring{end+1}='\sigma_{1}';

%------------------------
divider=11;
x_=logspace(floor(log10(min(FluctuationList))),ceil(log10(max(FluctuationList))),100);
h=plot(x_,x_./divider,'-k','linewidth',1);

x=FluctuationList(FluctuationList<0.05);
y=S0_mean(FluctuationList<0.05);

legendindex(end+1)=h;
legendstring{end+1}=sprintf('\\sigma_{est} (\\delta<<1) = \\delta/%i',divider);

%-----------------------------
% %fit curve S1 polynomial with high ranges only
y=S0_mean(FluctuationList>=10);
p1=round(nanmean(y)*1e4)/1e4;
x=log10(FluctuationList);
h=plot(10.^[min(x) max(x)],[p1 p1],':k','linewidth',1); 
legendindex(end+1)=h;
legendstring{end+1}=sprintf('\\sigma_{est} (\\delta>>1) = %.4f',p1);

%------------------------
legend(legendindex,legendstring)

box on
%grid on
xlabel('Fluctuation Level Percentage (\delta)')
ylabel('Average Standard Deviation Measured ( \sigma_m )')
set(gca,'XScale','Log');
set(gca,'YScale','Log');
set(gca,'Color',[1 1 1]);
set(gcf,'Color','white')
ylim([0 0.25])

clear y x x_ p1 h Yfit_ Yresid_ SSresid_ SStotal_ rsqp_ legendindex legendstring S0_mean S1_mean

%% SHOW PROBABILITY CURVE
figure(5)
hold on
prob_mean=[];
for k=1:length(FluctuationList)
    prob_mean(k)=nanmean(estimated_alpha(k,:));
end
plot(FluctuationList, prob_mean,'o-r','markersize',6,'linewidth',2);
ylim([0 1])
if length(FluctuationList)>1
    xlim([FluctuationList(end) FluctuationList(1)])
end
box on
xlabel('Fluctuation Level Percentage (\delta)')
ylabel('Average Probability of True Ion Composition ( P_\alpha  )')
set(gca,'XScale','Log');
set(gca,'Color',[1 1 1]);
set(gcf,'Color','white')

clear h1 h2 prob_mean prob_std k















                
%% CLEAR UNUSED VARIABLES
%clear all
clearvars -except save_filename
load (sprintf(save_filename))                


%% SHOW ALL VALUES OF DIFFERENT FLUCTUATIONS AND INVALID REGRESSION 

figure(6)
hold on    

for k=1:length(FluctuationList)
    FluctuationPercentage=FluctuationList(k);

    %% Load File
    Number_PlasmaParameters = length(PlasmaParameters);
    filename=sprintf('%ip_%.1ernd_%irep_%iparam_1segm0.10over.mat',Number_PlasmaParameters,FluctuationPercentage,Number_Random_Realizations,Number_InputParameters_Expected);
    filename
    if exist(filename, 'file') == 2
        load(filename);
    else
        warning ('Warning: file "%s" does not exist!!\n',filename)
        continue
    end
    clear filename

    %Reload Number of Stored Parameters
    Number_InputParameters=length(LM_results.fit_value);

    
    %% True Parameters Selection
    load('True_Rand_Values.mat')

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

    
    %% Filter Chi-square Minimum Limit

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
    
    %% Show Error vs True Parameter in different colors for different fluctuations
    
    X=LM_results.true_parameters(:,:,Selection_Parameter_index);
    Y=LM_results_valid.fit_parameters(:,:,Selection_Parameter_index);
    
    plot(X,Y,'.','markersize',4,'color',Fluctuation_Colors(k,:))
    box on
    grid on
    xlabel('True Ion Composition')
    ylabel('Fitted Ion Composition')
    ylim([0 1])
    set(gca,'Color',[1 1 1]);
    set(gcf,'Color','white')
    
end


%--------------------------------------------------------------------------
%Colorbar
cbh=colorbar();
liststr={};
for k=1:length(FluctuationList)
    if FluctuationList(k)<0.01
        liststr{k}=sprintf('%.3f',FluctuationList(k));
    elseif FluctuationList(k)<0.1
        liststr{k}=sprintf('%.2f',FluctuationList(k));
    elseif FluctuationList(k)<1
        liststr{k}=sprintf('%.1f',FluctuationList(k));
    else
        liststr{k}=sprintf('%i',FluctuationList(k));
    end
end
set(cbh,'ylim',[1 65], 'YTick',linspace(1,65, length(FluctuationList)),'YTickLabels',liststr,'YGrid','on')
ylabel(cbh,'Fluctuation Percentage \delta(%)','rotation',270)
%--------------------------------------------------------------------------
clear cbh liststr k S1 M1 S0 M0 X Y i minindex Xsize

%% Show Tendency previously calculated

hold on

X=[0:0.01:1];
a=-0.96873;
b=1.114298;
Y=a.*X+b;
plot(X,Y,'-k','linewidth',2)
Y=(a).*X+(b+0.18);
plot(X,Y,'--k','linewidth',1)
Y=(a).*X+(b-0.18);
plot(X,Y,'--k','linewidth',1)

clear X Y a b

%% Weighted Linear Regression 
figure(6)
hold on    

X=LM_results.true_parameters(:,:,Selection_Parameter_index);
X=sort(X(:,1));
X__=X(:,1);

M1=estimated_mean_1(1:end,:);
S1=estimated_sigma_1(1:end,:);
S1inv=1./S1;

M1=M1(9:end-4,:);
S1inv=S1inv(9:end-4,:);

S1inv(S1inv==Inf)=NaN;

%Weighed average of all M1 points at every X step
M1weighted=[];
for i=1:length(M1)
    M1weighted(i)=nansum(S1inv(:,i).*(M1(:,i)))./nansum(S1inv(:,i));
end

x=X__(~isnan(M1weighted));
y=x-M1weighted(~isnan(M1weighted))';

y=y(~isnan(x));
x=x(~isnan(x));

[InvalidFitResult,ErrorEst] = polyfit(x,y,1);
[YFit,Ydelta]=polyval(InvalidFitResult,X__,ErrorEst);

hold on
plot(X__,YFit,'-','markersize',5,'color','k','linewidth',2)        
plot(X__,YFit+2*Ydelta,'--','markersize',5,'color','k','linewidth',1)        
plot(X__,YFit-2*Ydelta,'--','markersize',5,'color','k','linewidth',1)  
ylim([0 1])
title(sprintf('Invalid Fitting: p_{invalid} = %.6f·p_{true}+%.6f (\\pm%.6f)',InvalidFitResult(1),InvalidFitResult(2),2*nanmean(Ydelta)))
InvalidFitResult

clear i Xsize X__ S1 S1inv M1  x y YFit Ydelta ErrorEst
clear M1weighted




















                
%% CLEAR UNUSED VARIABLES
%clear all
clearvars -except save_filename
load (sprintf(save_filename))                

%% SHOW PARAMETER ERROR VALUES FOR DIFFERENT FLUCTUATIONS (PERCENTAGE)

figure(7)
hold on    

newFluctuationList=FluctuationList(3:end);

Fluctuation_Colors=jet(length(newFluctuationList));

for k=1:length(newFluctuationList)
    FluctuationPercentage=newFluctuationList(k);

    %% Load File
    Number_PlasmaParameters = length(PlasmaParameters);
    filename=sprintf('%ip_%.1ernd_%irep_%iparam_1segm0.10over.mat',Number_PlasmaParameters,FluctuationPercentage,Number_Random_Realizations,Number_InputParameters_Expected);
    filename
    if exist(filename, 'file') == 2
        load(filename);
    else
        warning ('Warning: file "%s" does not exist!!\n',filename)
        continue
    end
    clear filename

    %Reload Number of Stored Parameters
    Number_InputParameters=length(LM_results.fit_value);

    
    %% True Parameters Selection
    load('True_Rand_Values.mat')

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

    clear True_Parameters_Vector True_Parameters ListOfPlasmaParameters 

    
    %% Filter Chi-square Minimum Limit

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
    
    %% Show Error vs True Parameter in different colors for different fluctuations
    for i=1:Number_PlasmaParameters-1
        subplot (1,Number_PlasmaParameters-1,i)
        hold on

        X=LM_results.err_parameters(:,:,Selection_Parameter_index);
        Y=100.*LM_results_valid.err_parameters(:,:,i)./LM_results_valid.true_parameters(:,:,i);

        plot(X,Y,'.','markersize',4,'color',Fluctuation_Colors(k,:))
        box on
        grid on
        xlabel('Error in Ion Composition')
        
        parIndex=find(ListOfPlasmaParametersIndex==i);
        if parIndex==1
            ylabel('Error in Electron Density (%)')
            ylim([-100 100])
        elseif parIndex==2
            ylabel('Error in Electron Temperature (%)')
            ylim([-200 100])
        elseif parIndex==3
            ylabel('Error in Ion Temperature (%)')
            ylim([-200 100])
        elseif parIndex==4
            ylabel('Error in Ion Drift Velocity (%)')
            ylim([-100 100])
        elseif parIndex==5
            ylabel('Error in Ion Composition (%)')
            ylim([-100 100])
        end
        
        set(gcf,'Color','white')
        set(gcf,'Position', [100, 100, 800, 600])

    end
    
end

%--------------------------------------------------------------------------
%Colorbar
for k=1:Number_PlasmaParameters-1
    subplot (1,Number_PlasmaParameters-1,k)

    cbh=colorbar();
    liststr={};
    for k=1:length(newFluctuationList)
        if newFluctuationList(k)<0.01
            liststr{k}=sprintf('%.3f',newFluctuationList(k));
        elseif newFluctuationList(k)<0.1
            liststr{k}=sprintf('%.2f',newFluctuationList(k));
        elseif newFluctuationList(k)<1
            liststr{k}=sprintf('%.1f',newFluctuationList(k));
        else
            liststr{k}=sprintf('%i',newFluctuationList(k));
        end
    end
    set(cbh,'ylim',[1 65], 'YTick',linspace(1,65, length(newFluctuationList)),'YTickLabels',liststr,'YGrid','on')
    ylabel(cbh,'Fluctuation Percentage \delta(%)','rotation',270)
end
%--------------------------------------------------------------------------
clear cbh liststr k S1 M1 S0 M0 X Y i minindex Fluctuation_Colors Xsize parIndex ListOfPlasmaParametersIndex









%% SHOW PARAMETER ERROR VALUES FOR DIFFERENT FLUCTUATIONS (PERCENTAGE) HISTOGRAM

figure(8)
hold on    

newFluctuationList=FluctuationList(3:end);

Fluctuation_Colors=jet(length(newFluctuationList));

for k=1:length(newFluctuationList)
    FluctuationPercentage=newFluctuationList(k);

    %% Load File
    Number_PlasmaParameters = length(PlasmaParameters);
    filename=sprintf('%ip_%.1ernd_%irep_%iparam_1segm0.10over.mat',Number_PlasmaParameters,FluctuationPercentage,Number_Random_Realizations,Number_InputParameters_Expected);
    filename
    if exist(filename, 'file') == 2
        load(filename);
    else
        warning ('Warning: file "%s" does not exist!!\n',filename)
        continue
    end
    clear filename

    %Reload Number of Stored Parameters
    Number_InputParameters=length(LM_results.fit_value);

    
    %% True Parameters Selection
    load('True_Rand_Values.mat')

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

    clear True_Parameters_Vector True_Parameters ListOfPlasmaParameters 

    
    %% Filter Chi-square Minimum Limit

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
    
    
    %% Show Error vs True Parameter in different colors for different fluctuations
    for i=1:Number_PlasmaParameters
        subplot (Number_PlasmaParameters,1,i)
        hold on
        num_points=1000;
        
        parIndex=find(ListOfPlasmaParametersIndex==i);
        if parIndex==1
            X=linspace(-100,100,num_points);
        elseif parIndex==2
            X=linspace(-200,100,num_points);
        elseif parIndex==3
            X=linspace(-200,100,num_points);
        elseif parIndex==4
            X=linspace(-100,100,num_points);
        elseif parIndex==5
            X=linspace(-100,100,num_points);
        end
        
        Y=reshape(100.*LM_results_valid.err_parameters(:,:,i)./LM_results_valid.true_parameters(:,:,i),[Number_InputParameters*Number_Random_Realizations 1]);
        Y=hist(Y,X);
        %Delete peaks at ending points (p=0 and p=1)
        Y(1:2)=0;
        Y(end-1:end)=0;
        Y=100.*Y./nanmax(Y);

        
        Y=[Y 0 0];
        Y(Y==0)=1e-4;
        X=[X max(X) min(X)];
        fill(X,Y,Fluctuation_Colors(k,:),'LineStyle','none')
        plot(X,Y,'-k','linewidth',1)
        
        box on
        grid on
        
        if parIndex==1
            xlabel('Error in Electron Density (%)')
        elseif parIndex==2
            xlabel('Error in Electron Temperature (%)')
        elseif parIndex==3
            xlabel('Error in Ion Temperature (%)')
        elseif parIndex==4
            xlabel('Error in Ion Drift Velocity (%)')
        elseif parIndex==5
            xlabel('Error in Ion Composition (%)')
        end
        
        ylabel('Histogram (%)')
        set(gca,'Color',[1 1 1]);
        set(gcf,'Color','white')

        set(gcf,'Position', [100, 100, 800, 600])
    end
    
end


clear cbh liststr k S1 M1 S0 M0 X Y i minindex Fluctuation_Colors Xsize num_points parIndex ListOfPlasmaParametersIndex





