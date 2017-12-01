% show_valid_and_correct_probabilities.m
%
% -------------------------------------------------------------------------
% Author: Miguel Martinez Ledesma (University of Chile)  
% Email: miguel.martinez@ing.uchile.cl
% Date: 24 November 2016
% -------------------------------------------------------------------------

%% Initialize
clear all


filenames={ 'Analyzed_Results_5p.mat'; ...
            'Analyzed_Results_4p_given_Vi.mat'; ...
            'Analyzed_Results_3p_given_NeVi.mat'; ...
            'Analyzed_Results_2p_given_NeViTe_Ti.mat'; ...
            'Analyzed_Results_2p_given_NeViTe.mat'};

colors=[];
colors(1,:)=[1 0 0];
colors(2,:)=[0 1 0];
colors(3,:)=[1 0.5 0];
colors(4,:)=[1 0 1];
colors(5,:)=[0 0 1];

symbols={};
symbols{1}='o-';
symbols{2}='+-';
symbols{3}='s-';
symbols{4}='p-';
symbols{5}='d-';

legendstring={};
legendstring{1}='5 param. (N_e, T_e, T_i, V_i, p)';
legendstring{2}='4 param. (N_e, T_e, T_i, p) given V_i';
legendstring{3}='3 param. (T_e, T_i, p) given N_e and V_i';
legendstring{4}='2 param. (T_i, p) given N_e, V_i, and T_e/T_i';
legendstring{5}='2 param. (T_i, p) given N_e, V_i, and T_e';

%% Plot
figure

for i=1:length(filenames)
    
    load(filenames{i})

    %%
    
    subplot(3,1,1)
    hold on
    plot(FluctuationList,100*Solved_all/(Number_InputParameters*Number_Random_Realizations),symbols{i},'color',colors(i,:),'markersize',6,'linewidth',1.5)
    grid on
    box on
    set(gca,'XScale','Log');
    xlabel('Fluctuation Percentage ( \delta(%) )')
    ylabel('P_{valid} (%)')
    legend(legendstring,'Location','Best','FontSize',9)
    ylim([40 100])

    subplot(3,1,2)
    hold on
    plot(FluctuationList,100*Correct_all./(Solved_all),symbols{i},'color',colors(i,:),'markersize',6,'linewidth',1.5)
    grid on
    box on
    set(gca,'XScale','Log');
    xlabel('Fluctuation Percentage ( \delta(%) )')
    ylabel('P_{correct} (%)')
    ylim([40 100])

    subplot(3,1,3)
    hold on
    plot(FluctuationList,100*Correct_all./(Number_InputParameters*Number_Random_Realizations),symbols{i},'color',colors(i,:),'markersize',6,'linewidth',1.5)
    grid on
    box on
    set(gca,'XScale','Log');
    xlabel('Fluctuation Percentage ( \delta(%) )')
    ylabel('P_{Valid & Correct} (%)')
    %legend(legendstring,'Location','Best')
    ylim([40 100])
end

set(gcf,'Position', [50, 0, 800, 1000])
