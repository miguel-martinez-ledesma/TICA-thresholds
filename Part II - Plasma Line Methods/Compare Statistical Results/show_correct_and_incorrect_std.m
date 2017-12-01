% show_correct_and_incorrect_std.m
%
% -------------------------------------------------------------------------
% Author: Miguel Martínez Ledesma (University of Chile)  
% Email: miguel.martinez@ing.uchile.cl
% Date: 24 November 2016
% -------------------------------------------------------------------------

%% Initialize
clear all


filenames={ 'Graphically_Analyzed_Results_5p.mat'; ...
            'Graphically_Analyzed_Results_4p_given_Vi.mat'; ...
            'Graphically_Analyzed_Results_3p_given_NeVi.mat'; ...
            'Graphically_Analyzed_Results_2p_given_NeViTe_Ti.mat'; ...
            'Graphically_Analyzed_Results_2p_given_NeViTe.mat'};


string={};
string{1}='5 param. (N_e, T_e, T_i, V_i, p)';
string{2}='4 param. (N_e, T_e, T_i, p) given V_i';
string{3}='3 param. (T_e, T_i, p) given N_e and V_i';
string{4}='2 param. (T_i, p) given N_e, V_i, and T_e/T_i';
string{5}='2 param. (T_i, p) given N_e, V_i, and T_e';

%% Plot
figure

for j=1:length(filenames)
    
    %%
	load(filenames{j})

    %%  
    subplot(3,2,j)
    hold on

    legendindex=[];
    legendstring={};

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
    h=plot(FluctuationList, S0_mean,'o-r','markersize',6,'linewidth',1.5);
    legendindex(end+1)=h;
    legendstring{end+1}='\sigma_{0}';

    %------------------------
    h=plot(FluctuationList(1:length(S1_mean)), S1_mean,'s-','color',[0.1 0.8 0.1] ,'markersize',6,'linewidth',1.5);
    legendindex(end+1)=h;
    legendstring{end+1}='\sigma_{1}';

    %------------------------
    switch j
        case 1
            divider=6;
        case 2 
            divider=6;
        case 3 
            divider=11;
        case 4
            divider=15;
        case 5
            divider=100;
        otherwise
            divider=15;
    end
    x_=logspace(floor(log10(min(FluctuationList))),ceil(log10(max(FluctuationList))),100);
    h=plot(x_,x_./divider,'--b','linewidth',1.5);

    x=FluctuationList(FluctuationList<0.05);
    y=S0_mean(FluctuationList<0.05);
    
    legendindex(end+1)=h;
    legendstring{end+1}=sprintf('\\sigma_{est} (\\delta<<1) = \\delta/%.0f',divider);

    %-----------------------------
    %fit curve S1 polynomial with high ranges only
    y=S0_mean(FluctuationList>=50);
    p1=round(nanmean(y)*1e4)/1e4;
    x=log10(FluctuationList);
    h=plot(10.^[min(x) max(x)],[p1 p1],':k','linewidth',1.5); 
    legendindex(end+1)=h;
    legendstring{end+1}=sprintf('\\sigma_{est} (\\delta>>1) = %.4f',p1);

    %------------------------
    legend(legendindex,legendstring,'Location','Best','FontSize',9)

    box on
    %grid on
    xlabel('Fluctuation Level Percentage ( \delta(%) )')
    ylabel('Standard Deviation ( \sigma )')
    set(gca,'XScale','Log');
    set(gca,'YScale','Log');
    set(gca,'Color',[1 1 1]);
    set(gcf,'Color','white')
    ylim([0 0.25])
    title(string{j})
    
    clear y x x_ p1 h Yfit_ Yresid_ SSresid_ SStotal_ rsqp_ S0_mean S1_mean legendstring legendindex


end


set(gcf,'Position', [50, 0, 1000, 1000])
