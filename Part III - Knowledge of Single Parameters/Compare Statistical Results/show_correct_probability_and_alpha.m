% show_correct_probability_and_alpha.m
%
% -------------------------------------------------------------------------
% Author: Miguel Martinez Ledesma (University of Chile)  
% Email: miguel.martinez@ing.uchile.cl
% Date: 24 November 2016
% -------------------------------------------------------------------------

%% Initialize
clear all


filenames_prob={ 'Graphically_Analyzed_Results_4p_given_Ne.mat'; ...
            'Graphically_Analyzed_Results_4p_given_Te_Ti.mat'; ...
            'Graphically_Analyzed_Results_4p_given_Te.mat'; ...
            'Graphically_Analyzed_Results_4p_given_Ti.mat';...
            'Graphically_Analyzed_Results_4p_given_Vi.mat';...
            'Graphically_Analyzed_Results_5p.mat'};
        
        
filenames_stat={ 'Analyzed_Results_4p_given_Ne.mat'; ...
            'Analyzed_Results_4p_given_Te_Ti.mat'; ...
            'Analyzed_Results_4p_given_Te.mat'; ...
            'Analyzed_Results_4p_given_Ti.mat';...
            'Analyzed_Results_4p_given_Vi.mat';...
            'Analyzed_Results_5p'};

string={};
string{1}='4 param. (T_e, T_i, V_i, p) given N_e';
string{2}='4 param. (N_e, T_i, V_i, p) given T_e/T_i';
string{3}='4 param. (N_e, T_e, V_i, p) given T_e';
string{4}='4 param. (N_e, T_i, V_i, p) given T_i';
string{5}='4 param. (N_e, T_e, T_i, p) given V_i';
string{6}='5 param. (N_e, T_e, T_i, V_i, p)';

colors=[];
colors(1,:)=[1 0.5 0.5];
colors(2,:)=[0 1 0];
colors(3,:)=[1 0 1];
colors(4,:)=[0 0.5 0.5];
colors(5,:)=[0 0 1];
colors(6,:)=[1 0 0];
colors(7,:)=[0.5 0.5 0];

symbols_stat={};
symbols_stat{1}='o-';
symbols_stat{2}='d-';
symbols_stat{3}='s-';
symbols_stat{4}='p-';
symbols_stat{5}='^-';
symbols_stat{6}='+-';
symbols_stat{7}='.-';

symbols_prob={};
symbols_prob{1}='--';
symbols_prob{2}='--';
symbols_prob{3}='--';
symbols_prob{4}='--';
symbols_prob{5}='--';
symbols_prob{6}='--';
symbols_prob{7}='--';


linewidth=[];
linewidth(1)=1.5;
linewidth(2)=1.5;
linewidth(3)=1.5;
linewidth(4)=1.5;
linewidth(5)=1.5;
linewidth(6)=1.5;
linewidth(7)=1.5;

%% Plot
figure    

for j=1:length(filenames_stat)
    
    subplot(3,round(length(filenames_stat)/3),j)
    
    
    %%
    load(filenames_stat{j})

    %%
    hold on
    plot(FluctuationList,Correct_all./(Solved_all),symbols_stat{j},'color',colors(j,:),'markersize',6,'linewidth',linewidth(j));
    grid on
    box on
    set(gca,'XScale','Log');
    xlabel('Fluctuation Percentage ( \delta(%) )')
    ylabel('Probability of Correct Ion Composition')
    set(gca,'Color',[1 1 1]);
    set(gcf,'Color','white')
    title(string{j})
    ylim([0.30 1])
    
    %%
    load(filenames_prob{j})


    %%
    hold on
    prob_mean=[];
    for k=1:length(FluctuationList)
        prob_mean(k)=nanmean(estimated_alpha(k,:));
    end
    plot(FluctuationList, prob_mean,symbols_prob{j},'color',colors(j,:),'markersize',6,'linewidth',linewidth(j));
    ylim([0.3 1])
    xlim([FluctuationList(end) FluctuationList(1)])

    box on
    grid on
    xlabel('Fluctuation Level Percentage (\delta)')
    ylabel('Probability of Correct Ion Composition')
    set(gca,'XScale','Log');
    set(gca,'Color',[1 1 1]);
    set(gcf,'Color','white')
    
    clear h1 h2 prob_mean prob_std k

    legend({'P_{correct}'; 'Average Weight \alpha'},'Location','Best')
    

end


set(gcf,'Position', [10, 10, 1000, 900])
    