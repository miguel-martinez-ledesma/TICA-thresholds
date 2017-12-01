% show_valid_and_correct_probabilities.m
%
% -------------------------------------------------------------------------
% Author: Miguel Martinez Ledesma (University of Chile)  
% Email: miguel.martinez@ing.uchile.cl
% Date: 24 November 2016
% -------------------------------------------------------------------------

%% Initialize
clear all


filenames={ 'Analyzed_Results_4p_given_Ne.mat'; ...
            'Analyzed_Results_4p_given_Te_Ti.mat'; ...
            'Analyzed_Results_4p_given_Te.mat'; ...
            'Analyzed_Results_4p_given_Ti.mat';...
            'Analyzed_Results_4p_given_Vi.mat';...
            'Analyzed_Results_5p'};

colors=[];
colors(1,:)=[1 0.5 0.5];
colors(2,:)=[0 1 0];
colors(3,:)=[1 0 1];
colors(4,:)=[0 0.5 0.5];
colors(5,:)=[0 0 1];
colors(6,:)=[1 0 0];
colors(7,:)=[0.5 0.5 0];

symbols={};
symbols{1}='o-';
symbols{2}='d-';
symbols{3}='s-';
symbols{4}='p-';
symbols{5}='^-';
symbols{6}='+--';
symbols{7}='.--';

linewidth=[];
linewidth(1)=1.5;
linewidth(2)=1.5;
linewidth(3)=1.5;
linewidth(4)=1.5;
linewidth(5)=1.5;
linewidth(6)=1.5;
linewidth(7)=1.5;

legendstring={};
legendstring{1}='4 param. (T_e, T_i, V_i, p) given N_e';
legendstring{2}='4 param. (N_e, T_i, V_i, p) given T_e/T_i';
legendstring{3}='4 param. (N_e, T_i, V_i, p) given T_e';
legendstring{4}='4 param. (N_e, T_e, V_i, p) given T_i';
legendstring{5}='4 param. (N_e, T_e, T_i, p) given V_i';
legendstring{6}='5 param. (N_e, T_e, T_i, V_i, p)';

%% Plot
figure

for i=1:length(filenames)
    
    load(filenames{i})

    %%
    
    subplot(3,1,1)
    hold on
    plot(FluctuationList,100*Solved_all/(Number_InputParameters*Number_Random_Realizations),symbols{i},'color',colors(i,:),'markersize',6,'linewidth',linewidth(i))
    grid on
    box on
    set(gca,'XScale','Log');
    xlabel('Fluctuation Percentage ( \delta(%) )')
    ylabel('P_{valid} (%)')
    legend(legendstring,'Location','Best','FontSize',9)
    ylim([40 100])

    subplot(3,1,2)
    hold on
    plot(FluctuationList,100*Correct_all./(Solved_all),symbols{i},'color',colors(i,:),'markersize',6,'linewidth',linewidth(i))
    grid on
    box on
    set(gca,'XScale','Log');
    xlabel('Fluctuation Percentage ( \delta(%) )')
    ylabel('P_{correct} (%)')
    %legend(legendstring,'Location','Best')
    ylim([40 100])

    subplot(3,1,3)
    hold on
    plot(FluctuationList,100*Correct_all./(Number_InputParameters*Number_Random_Realizations),symbols{i},'color',colors(i,:),'markersize',6,'linewidth',linewidth((i)))
    grid on
    box on
    set(gca,'XScale','Log');
    xlabel('Fluctuation Percentage ( \delta(%) )')
    ylabel('P_{valid & correct} (%)')
    %legend(legendstring,'Location','Best')
    ylim([40 100])
end

set(gcf,'Position', [50, 0, 800, 1000])
