% Description: this script makes the figures used in the final main body of
% the manuscript

% Update to your local path
dir = '/Users/keaniw/Documents/Research/Kenya Project/Project Code/Flexible_Ops_Plan_Design/';
addpath(genpath(dir))

%% Figure 3: Mean simulated total costs by infrastructure alternative

% Description: Mean simulated total costs by infrastructure alternative. 
% From the total 10,000 Monte Carlo simulations, we plot the mean simulated 
% total costs, mean initial dam costs, and mean expansion costs for each 
% flexible infrastructure alternative subset by simulations that transition
% to representative dry (66-76 mm/mo), moderate (77-86 mm/mo), and wet 
% (87-97 mm/mo) final climate states.

% =========================================================================
% (1) LOAD THE DATA FILES TO PLOT

% specify file naming convensions
discounts = 0;
cprimes = 1.25e-6;

cd('/Users/keaniw/Documents/Research/Kenya Project/Project Code/Flexible_Ops_Plan_Design/')

% specify file path for data
folder = 'Results/Results_SDP_expansion';


for d=1:length(discounts)
    disc = string(discounts(d));
    for c = 1:length(cprimes)
        cp = cprimes(c);
        c_prime = regexprep(strrep(string(cp), '.', ''), {'-0'}, {''});
        
        % load adaptive operations files:
        load(strcat(folder,'/BestFlex_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestStatic_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestPlan_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestFlexStaticPlan_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))

        bestAct_adapt = bestAct;
        V_adapt = V;
        Vs_adapt = Vs;
        Vd_adapt = Vd;
        Vs_static_adapt = allVs_static;
        Vs_flex_adapt = allVs_flex;
        Vs_plan_adapt = allVs_plan;
        Vd_static_adapt = allVd_static;
        Vd_flex_adapt = allVd_flex;
        Vd_plan_adapt = allVd_plan;
        X_adapt = X;
        C_adapt = C_state;
        action_adapt = action;
        totalCostTime_adapt = totalCostTime;
        damCostTime_adapt = damCostTime;
        P_state_adapt = P_state;
        bestVal_flex_adapt = bestVal_flex;
        bestVal_static_adapt = bestVal_static;
        bestVal_plan_adapt = bestVal_plan;
        if bestAct(3) + bestAct(4)*bestAct(5) > 150
            bestAct_adapt(4) = (150 - bestAct(2))/bestAct(5);
        end
        if bestAct(8) + bestAct(9)*bestAct(10) > 150
            bestAct_adapt(9) = (150 - bestAct(8))/bestAct(10);
        end
        
        % load non-adaptive operations files:
        load(strcat(folder,'/BestFlex_nonadaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestStatic_nonadaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestPlan_nonadaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestFlexStaticPlan_nonadaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        
        bestAct_nonadapt = bestAct;
        V_nonadapt = V;
        Vs_nonadapt = Vs;
        Vd_nonadapt = Vd;
        Vs_static_nonadapt = allVs_static;
        Vs_flex_nonadapt = allVs_flex;
        Vs_plan_nonadapt = allVs_plan;
        Vd_static_nonadapt = allVd_static;
        Vd_flex_nonadapt = allVd_flex;
        Vd_plan_nonadapt = allVd_plan;
        X_nonadapt = X;
        C_nonadapt = C_state;
        action_nonadapt = action;
        totalCostTime_nonadapt = totalCostTime;
        damCostTime_nonadapt = damCostTime;
        P_state_nonadapt = P_state;
        bestVal_flex_nonadapt = bestVal_flex;
        bestVal_static_nonadapt = bestVal_static;
        bestVal_plan_nonadapt = bestVal_plan;
        if bestAct(3) + bestAct(4)*bestAct(5) > 150
            bestAct_nonadapt(4) = (150 - bestAct(2))/bestAct(5);
        end
        if bestAct(8) + bestAct(9)*bestAct(10) > 150
            bestAct_nonadapt(9) = (150 - bestAct(8))/bestAct(10);
        end
        
    end
      
    % =====================================================================
    % CREATE INFRASTRUCTURE COST LOOK UP STRUCTURES
    % Note: for discounted infrastructure and shortage costs, use: totalCostTime and damCostTime
    
    % calculate non-discounted infrastructure costs using storage2damcost:
    for s = 1:2 % for non-adaptive and adaptive operations
        if s == 1
            x = bestAct_nonadapt;
        else
            x = bestAct_adapt;
        end
        
        % load optimal infrastructure information
        optParam.staticCap = x(2); % static dam size [MCM]
        optParam.smallFlexCap = x(3); % unexpanded flexible design dam size [MCM]
        optParam.numFlex = x(4);  % number of possible expansion capacities [#]
        optParam.flexIncr = x(5); % increment of flexible expansion capacities [MCM]
        costParam.PercFlex = x(6); % Initial upfront capital cost increase (0.075)
        costParam.PercFlexExp = x(7); % Expansion cost of flexible dam  (0.15)
        optParam.smallPlanCap = x(8); % unexpanded flexible planning dam size [MCM]
        optParam.numPlan = x(9); % MCM
        optParam.planIncr = x(10);
        costParam.PercPlan = x(11); % initial upfront capital cost increase (0);
        costParam.PercPlanExp = x(12); % expansion cost of flexibly planned dam (0.5)
        costParam.discountrate = x(15);
        
        s_C = 1:3+optParam.numFlex+optParam.numPlan;
        M_C = length(s_C);
        
        storage = zeros(1, M_C);
        storage(1) = optParam.staticCap;
        storage(2) = optParam.smallFlexCap;
        storage(3) = optParam.smallPlanCap;
        storage(4:3+optParam.numFlex) = min(storage(2) + (1:optParam.numFlex)*optParam.flexIncr, 150);
        storage(4+optParam.numFlex:end) = min(storage(3) + (1:optParam.numPlan)*optParam.planIncr, 150);
        
        % Actions: Choose dam option in time period 1; 
        % expand dam in future time periods
        a_exp = 0:3+optParam.numFlex+optParam.numPlan;
        
        % Define infrastructure costs
        infra_cost = zeros(1,length(a_exp));
        infra_cost(2) = storage2damcost(storage(1),0); % cost of static dam
        for i = 1:optParam.numFlex % cost of flexible design dam
            [infra_cost(3), infra_cost(i+4)] = storage2damcost(storage(2), ...
                storage(i+3),costParam.PercFlex, costParam.PercFlexExp); % cost of flexible design exp to option X
        end
        for i = 1:optParam.numPlan % cost of flexible plan dam
            [infra_cost(4), infra_cost(i+4+optParam.numFlex)] = storage2damcost(storage(3), ...
                storage(i+3+optParam.numFlex),costParam.PercPlan, costParam.PercPlanExp); % cost of flexible planning exp to option X
        end
        
        if s == 1
            infra_cost_nonadaptive = struct();
            infra_cost_nonadaptive.static = [storage(1); infra_cost(2)/1E6];
            infra_cost_nonadaptive.flex = [[storage(2), storage(4:3+optParam.numFlex)]; [infra_cost(3), infra_cost(5:optParam.numFlex+4)]./1E6];
            infra_cost_nonadaptive.plan = [[storage(3), storage(4+optParam.numFlex:end)]; [infra_cost(4), infra_cost(5+optParam.numFlex:end)]./1E6];
            infra_cost_nonadaptive_lookup = infra_cost;
        else
            infra_cost_adaptive = struct();
            infra_cost_adaptive.static = [storage(1); infra_cost(2)/1E6];
            infra_cost_adaptive.flex = [[storage(2), storage(4:3+optParam.numFlex)]; [infra_cost(3), infra_cost(5:optParam.numFlex+4)]./1E6];
            infra_cost_adaptive.plan = [[storage(3), storage(4+optParam.numFlex:end)]; [infra_cost(4), infra_cost(5+optParam.numFlex:end)]./1E6];
            infra_cost_adaptive_lookup = infra_cost;
        end
    end

% =========================================================================

    % create bar plot of expected dam and shortage costs:
    f=figure('units','centimeters','position',[0,0,19,9]);
    t = tiledlayout(3,1);
    t.TileSpacing = 'tight';
    
    P_ranges = {[66:1:76];[77:1:86];[87:97]}; % reference dry, moderate, and wet climates
    
     s_C_adapt = [bestAct_adapt(2) bestAct_adapt(3), bestAct_adapt(8) (bestAct_adapt(3)+[1:bestAct_adapt(4)]*bestAct_adapt(5)),...
        (bestAct_adapt(8)+[1:bestAct_adapt(9)]*bestAct_adapt(10))];
    
    s_C_nonadapt = [bestAct_nonadapt(2) bestAct_nonadapt(3), bestAct_nonadapt(8) (bestAct_nonadapt(3)+[1:bestAct_nonadapt(4)]*bestAct_nonadapt(5)),...
        (bestAct_nonadapt(8)+[1:bestAct_nonadapt(9)]*bestAct_adapt(10))];
    
    % forward simulations of shortage and infrastructure costs
    totalCost_adapt = squeeze(totalCostTime_adapt(:,:,1:3))/1E6;% 1 is flex, 2 is static
    totalCost_nonadapt = squeeze(totalCostTime_nonadapt(:,:,1:3))/1E6;% 1 is flex, 2 is static
    damCost_adapt = squeeze(damCostTime_adapt(:,:,1:3))/1E6;% 1 is flex, 2 is static
    damCost_nonadapt = squeeze(damCostTime_nonadapt(:,:,1:3))/1E6;% 1 is flex, 2 is static
    
    for p = 1:length(P_ranges)
        % re-initialize mean cost arrays:
        totalCost_P = NaN(6,1);
        damCost_P = NaN(6,1);
        expCost_P = NaN(6,1);
        shortageCost_P = NaN(6,1);
        
        % calculate the average shortage, dam, and expansion costs for each
        % climate state:
        ind_P = ismember(P_state_adapt(:,5),P_ranges{p});
        
        % find climate subset from data for adaptive and non-adaptive
        % operations
        totalCost_P([2,1,3]) = squeeze(mean(sum(totalCost_nonadapt(ind_P,:,:),2)));
        damCost_P([2,1,3]) = squeeze(mean(damCost_nonadapt(ind_P,1,:)));
        expCost_P([2,1,3]) = squeeze(mean(sum(damCost_nonadapt(ind_P,2:5,:),2)));
        shortageCost_P([2,1,3]) = squeeze(mean(sum(totalCost_nonadapt(ind_P,2:5,:)- damCost_nonadapt(ind_P,2:5,:),2)));    
        
        totalCost_P([5,4,6]) = squeeze(mean(sum(totalCost_adapt(ind_P,:,:),2)));
        damCost_P([5,4,6]) = squeeze(mean(damCost_adapt(ind_P,1,:)));
        expCost_P([5,4,6]) = squeeze(mean(sum(damCost_adapt(ind_P,2:5,:),2)));
        shortageCost_P([5,4,6]) = squeeze(mean(sum(totalCost_adapt(ind_P,2:5,:)- damCost_adapt(ind_P,2:5,:),2)));    

        costs = [damCost_P, expCost_P, shortageCost_P];
        
        % bar plot of average costs from simulations:
        nexttile
        c = [[230,230,230];[110, 110, 110];[177,89,89]; [114,114,114];[90,90,90];[44,44,44]]/255;
        b = bar([1,2,3,4.2,5.2,6.2],costs,'stacked','FaceColor',"flat");
        
        if p == 1
            title({strcat("Dry Final Climates")},'FontSize',13)
            set(gca, 'XTick',[],'FontSize',0.1)
            h=gca; 
            h.XAxis.TickLength = [0 0];
            set(gca,'XTickLabel',[],'FontSize',0.1)
            leg = legend('Initial Dam Cost','Dam Expansion Cost','Shortage Cost','FontSize',6);
            set(leg,'Box','off')
        else
            title({strcat("Wet Final Climates")},'FontSize',13)
            if p==2
                title({strcat("Moderate Final Climates")},'FontSize',13)
                ylabel('Mean Simulated Cost (M$)','FontWeight','bold')
                set(gca, 'XTick',[],'FontSize',0.1)
                h=gca; 
                h.XAxis.TickLength = [0 0];
                set(gca,'XTickLabel',[],'FontSize',0.1)
            end        
        end

        for k = 1:size(b,2)
            b(k).CData = c(k,:);
        end
        
        ax = gca;
        ax.LineWidth = 0.5;

        xlim([0.5,6.7])
        yl = ylim;
        ylim([0, yl(2)+0.2*yl(2)])

        if p == 3

            set(gca, 'XTick', [1,2,3,4.2,5.2,6.2])
    
            set(gca,'XTickLabel',{strcat('Static Design'),...
                strcat('Flexible Design'),...
                strcat('Flexible Planning'),...
                strcat('Static Design'),...
                strcat('Flexible Design'),...
                strcat('Flexible Planning')})
            
            text(1.5, -35, "Static Operations",'FontSize',15,'FontWeight','bold')
            text(4.45, -35, "Flexible Operations",'FontSize',15,'FontWeight','bold')
                  
        end
    end
    
end

% save figure as .pdf and .eps and update fontsize
set(findall(gcf,'-property','FontSize'),'FontSize',8)
set(leg,'FontSize',6.5,'Position', leg.Position + [0.04 0.05 0 0])
cd(strcat(dir, 'Plots'))
exportgraphics(gcf, 'Figures/mainbody/fig3_barplot.pdf', 'ContentType', 'vector'); 
saveas(gcf, 'Figures/mainbody/fig3_barplot.eps')
saveas(gcf, 'Figures/mainbody/fig3_barplot.jpg')

%% Figure 4: Line plot of shortage cost vs. precipitation state

% Description: Plot of shortage costs vs. mean precipitation state for a 
% range of dam total installed capacities under flexible and static 
% operating policies.

% =========================================================================
% (1) LOAD THE DATA FILES TO PLOT

% We will subset results for target precip range of 66 to 97 mm/month
s_P_abs = 49:1:119; % expanded state space [mm/month]
storages = [50 70 90 110 130 150]; % storages to plot in MCM
cprime = 1.25E-6;

% set the shortagecost folder path
folder = 'Results/Results_SDP_reservoir_ops/';
scost_adapt = zeros(length(s_P_abs),length(storages));
scost_nonadapt = zeros(length(s_P_abs),length(storages));

% make data matrix for plots
for s=1:length(storages)

    s_now = storages(s);  % current dam capacity to load

    % load shortage costs for flexible operations
    load(strcat(folder,'sdp_adaptive_shortage_cost_s',string(s_now),'.mat'));    
    scost_adapt(:,s) = shortageCost(1,:).*cprime/1E6;
    
    % load shortage cost for static operaitons
    load(strcat(folder,'sdp_nonadaptive_shortage_cost_s',string(s_now),'.mat'));
    scost_nonadapt(:,s) = shortageCost(1,:).*cprime/1E6;
end

% =========================================================================
% (2) LINEPLOT: SHORTAGE COST VS. PRECIPITATION STATE

% initialize the figure graphics
f=figure('units','centimeters','position',[0,0,19,8]);
t = tiledlayout(1,2,'Padding','compact');
colors = {[0.6350 0.0780 0.1840] [0.8500 0.3250 0.0980] [0.9290 0.6940 0.1250]...
        [0.4660 0.6740 0.1880] [0 0.4470 0.7410] [0.4940 0.1840 0.5560]};


% make line plot for flexible and static operations
for i=1:2
    if i==1 % static operations
        p = plot(s_P_abs,scost_nonadapt,'LineWidth',1.2);
    else % flexible operations
        p = plot(s_P_abs,scost_adapt,'--','LineWidth',1.2);
    end
    
    % format graphics
    for j=1:length(storages)
        p(j).Color = colors{j};
    end
    xlabel('Precipitation State (mm/mo)','FontWeight','bold')
    ylabel('Shortage Cost (M$)','FontWeight','bold')
    xlim([66,80]);
    hold on
end

% add title and format legends
title('Shortage Cost vs. Precipitation State','FontWeight','bold')
leg1 = legend(strcat(string(storages), " MCM"),'Location','NorthEast', ...
    'AutoUpdate','off');
l1 = plot([-1, -2],[-1,-2],'color', [17 17 17]/255,'linestyle','-');
l2 = plot([-1, -2],[-1,-2],'color', [17 17 17]/255,'linestyle','--');
title(leg1,'Dam Capacity')
leg1.Title.Visible = 'on';
set(leg1,'Box','off')
ah1=axes('position',get(gca,'position'),'visible','off');
leg2 = legend(ah1, [l1 l2], "Static", "Flexible",'Location','SouthEast');
title(leg2,'Operations')
leg2.Title.Visible = 'on';
set(leg2,'Box','off','Position', leg2.Position + [0 0.2 0 0]);

% save figure as .pdf and .eps and update fontsize
set(findall(gcf,'-property','FontSize'),'FontSize',8)
cd(strcat(dir, 'Plots'))
exportgraphics(gcf, 'Figures/mainbody/fig4_shortagecosts.pdf', 'ContentType', 'vector'); 
saveas(gcf, 'Figures/mainbody/fig4_shortagecosts.eps')
saveas(gcf, 'Figures/mainbody/fig4_shortagecosts.jpg')

%% Figure 5: Box plot of cost savings for flex design and flex planning

% Description: Box plot of the precipitation state by time period across 
% simulations that yield the top 10% of total costs savings for flexible 
% planning and flexible design, given flexible operations 

% =========================================================================
% (1) LOAD THE DATA FILES TO PLOT

cd('/Users/keaniw/Documents/Research/Kenya Project/Project Code/Flexible_Ops_Plan_Design/')

% specify file pathc for data
folder = 'Results/Results_SDP_expansion';

% specify file naming convensions
disc = string(0);
cp = 1.25e-6;
c_prime = regexprep(strrep(string(cp), '.', ''), {'-0'}, {''});
         
% load flexible operations files
load(strcat(folder,'/BestFlexStaticPlan_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))

% specify variable names
totalCostTime_adapt = totalCostTime;
P_state_adapt = P_state;

% set different decade label parameters for use in plotting
decade_short = {'2001-20', '2021-40', '2041-60', '2061-80', '2081-00'};

% =========================================================================
% (2) PREPROCESS THE LOADED DATA

% make a generalized vector of 10,000 indices
allIdx = 1:length(totalCostTime_adapt);

% consider only flexible operations indices for flex planning and design
TotalCostDiff_adapt = (sum(totalCostTime_adapt(:,:,1),2) ...
    - sum(totalCostTime_adapt(:,:,3),2))/1E6; % cost flex design - flex plan

% find the row indices of top 10 percent of forward simulations where 
% flex design performs better than flexible planning (and vice-versa)
Idx_bestPlan_adapt = allIdx(TotalCostDiff_adapt >= ...
    prctile(TotalCostDiff_adapt(TotalCostDiff_adapt>0), 90));
Idx_bestFlex_adapt = allIdx(TotalCostDiff_adapt <= ...
    prctile(TotalCostDiff_adapt(TotalCostDiff_adapt<0), 10));

% =========================================================================
% (3) BOXPLOT: CLIMATE STATES OF  10% BEST SIMULATIONS FOR COST SAVINGS

% initialize the figure graphics
figure('units','centimeters','position',[0,0,19,6])
facecolors = [[153,204,204]; [187,221,221]; [153, 153, 204]; [204, 204, 255];...
    [255, 102, 102]; [255, 153, 153]]/255;

% make boxplot for simulations where flex design better than flex planning
boxplot([P_state_adapt(Idx_bestFlex_adapt,:)], 'Widths',0.2,'OutlierSize',5, ...
    'Symbol','.','BoxStyle','filled','positions', [(1:5)-0.2])

% format the graphics for boxplots where flex design better than flex plan
whisks = findobj(gca,'Tag','Whisker');
outs = findobj(gca, 'type', 'line','Tag', 'Outliers');
meds = findobj(gca, 'type', 'line', 'Tag', 'Median');
set(meds(1:5),'Color','k');
set(whisks(1:5),'Color',facecolors(4,:));
set(outs(1:5),'MarkerEdgeColor',facecolors(4,:));
a = findobj(gca,'Tag','Box');
set(a(1:5),'Color',facecolors(4,:),'LineWidth',20);

% make boxplot for simulations where flex planning better than flex design
hold on
boxplot([P_state_adapt(Idx_bestPlan_adapt,:)], 'Widths',0.2,'OutlierSize' ...
    ,5,'Symbol','.','BoxStyle','filled','positions', [(1:5)+0.2])

% format the graphics for boxplots where flex plan better than flex design
whisks = findobj(gca,'Tag','Whisker');
outs = findobj(gca, 'type', 'line','Tag', 'Outliers');
meds = findobj(gca, 'type', 'line', 'Tag', 'Median');
set(meds(1:5),'Color','k');
set(whisks(1:5),'Color',facecolors(6,:));
set(outs(1:5),'MarkerEdgeColor',facecolors(6,:));
a = findobj(gca,'Tag','Box');
set(a(1:5),'Color',facecolors(6,:),'LineWidth',20);

% plot mean value of cost savings for flex design and flex planning as dots
hold on
p1=plot((1:5)-0.2,mean(P_state_adapt(Idx_bestFlex_adapt,:)), 'k.');
plot((1:5)+0.2,mean(P_state_adapt(Idx_bestPlan_adapt,:)), 'k.')

% add figure legend
leg1 = legend([a(1) a(6)],{'Flexible Planning','Flexible Design'}, ...
    'Location','northwest');
set(leg1,'Box','off')

% format axes ticks, labels, and figure title
xticks(1:5)
xticklabels(decade_short);
ax=gca;
xlabel('Time Period','FontWeight','bold')
ylim([65,100])
set(gca, 'YGrid', 'off', 'YMinorGrid','off', 'XGrid', 'off');
ylabel('Precipitation State (mm/mo)','FontWeight','bold')
title({['Simulations with Greatest Cost Savings for Flexible ' ...
    'Planning vs. Flexible Design']},'Fontweight','bold')

% save figure as .pdf and .eps and update fontsize
set(findall(gcf,'-property','FontSize'),'FontSize',8)
cd(strcat(dir, 'Plots'))
exportgraphics(gcf, 'Figures/mainbody/fig5_boxplot.pdf', 'ContentType', 'vector'); 
saveas(gcf, 'Figures/mainbody/fig5_boxplot.eps')
saveas(gcf, 'Figures/mainbody/fig5_boxplot.jpg')


%% Figure 6: Mean simulated total costs for initially undersized dam

% Description: Mean simulated total costs by infrastructure alternative. 
% From the total 10,000 Monte Carlo simulations, we plot the mean simulated 
% total costs, mean initial dam costs, and mean expansion costs for each 
% flexible infrastructure alternative subset by simulations that transition
% to representative dry (66-76 mm/mo), moderate (77-86 mm/mo), and wet 
% (87-97 mm/mo) final climate states.
%
% Here, we assume reservoirs are initially sized to be 80 MCM, which is
% less than their previously found optimal capacities.

% =========================================================================
% (1) LOAD THE DATA FILES TO PLOT

% specify file naming convensions
discounts = 0;
cprimes = 1.25e-6;

cd('/Users/keaniw/Documents/Research/Kenya Project/Project Code/Flexible_Ops_Plan_Design/')

% specify file pathc for data
folder = 'Results/Results_SDP_expansion';


for d=1:length(discounts)
    disc = string(discounts(d));
    for c = 1:length(cprimes)
        cp = cprimes(c);
        c_prime = regexprep(strrep(string(cp), '.', ''), {'-0'}, {''});
        
        % load adaptive operations files:
        load(strcat(folder,'/BestFlex_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr120.mat'))
        load(strcat(folder,'/BestStatic_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr120.mat'))
        load(strcat(folder,'/BestPlan_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr120.mat'))
        load(strcat(folder,'/BestFlexStaticPlan_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr120.mat'))

        bestAct_adapt = bestAct;
        V_adapt = V;
        Vs_adapt = Vs;
        Vd_adapt = Vd;
        Vs_static_adapt = allVs_static;
        Vs_flex_adapt = allVs_flex;
        Vs_plan_adapt = allVs_plan;
        Vd_static_adapt = allVd_static;
        Vd_flex_adapt = allVd_flex;
        Vd_plan_adapt = allVd_plan;
        X_adapt = X;
        C_adapt = C_state;
        action_adapt = action;
        totalCostTime_adapt = totalCostTime;
        damCostTime_adapt = damCostTime;
        P_state_adapt = P_state;
        bestVal_flex_adapt = bestVal_flex;
        bestVal_static_adapt = bestVal_static;
        bestVal_plan_adapt = bestVal_plan;
        if bestAct(3) + bestAct(4)*bestAct(5) > 150
            bestAct_adapt(4) = (150 - bestAct(2))/bestAct(5);
        end
        if bestAct(8) + bestAct(9)*bestAct(10) > 150
            bestAct_adapt(9) = (150 - bestAct(8))/bestAct(10);
        end
        
        % load non-adaptive operations files:
        load(strcat(folder,'/BestFlex_nonadaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr120.mat'))
        load(strcat(folder,'/BestStatic_nonadaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr120.mat'))
        load(strcat(folder,'/BestPlan_nonadaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr120.mat'))
        load(strcat(folder,'/BestFlexStaticPlan_nonadaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr120.mat'))

        bestAct_nonadapt = bestAct;
        V_nonadapt = V;
        Vs_nonadapt = Vs;
        Vd_nonadapt = Vd;
        Vs_static_nonadapt = allVs_static;
        Vs_flex_nonadapt = allVs_flex;
        Vs_plan_nonadapt = allVs_plan;
        Vd_static_nonadapt = allVd_static;
        Vd_flex_nonadapt = allVd_flex;
        Vd_plan_nonadapt = allVd_plan;
        X_nonadapt = X;
        C_nonadapt = C_state;
        action_nonadapt = action;
        totalCostTime_nonadapt = totalCostTime;
        damCostTime_nonadapt = damCostTime;
        P_state_nonadapt = P_state;
        bestVal_flex_nonadapt = bestVal_flex;
        bestVal_static_nonadapt = bestVal_static;
        bestVal_plan_nonadapt = bestVal_plan;
        if bestAct(3) + bestAct(4)*bestAct(5) > 150
            bestAct_nonadapt(4) = (150 - bestAct(2))/bestAct(5);
        end
        if bestAct(8) + bestAct(9)*bestAct(10) > 150
            bestAct_nonadapt(9) = (150 - bestAct(8))/bestAct(10);
        end
        
    end
      
    % =====================================================================
    % CREATE INFRASTRUCTURE COST LOOK UP STRUCTURES
    % Note: for discounted infrastructure and shortage costs, use totalCostTime and damCostTime
    
    % calculate non-discounted infrastructure costs using storage2damcost:
    for s = 1:2 % for non-adaptive and adaptive operations
        if s == 1
            x = bestAct_nonadapt;
        else
            x = bestAct_adapt;
        end
        
        % load optimal infrastructure information
        optParam.staticCap = x(2); % static dam size [MCM]
        optParam.smallFlexCap = x(3); % unexpanded flexible design dam size [MCM]
        optParam.numFlex = x(4);  % number of possible expansion capacities [#]
        optParam.flexIncr = x(5); % increment of flexible expansion capacities [MCM]
        costParam.PercFlex = x(6); % Initial upfront capital cost increase (0.075)
        costParam.PercFlexExp = x(7); % Expansion cost of flexible dam  (0.15)
        optParam.smallPlanCap = x(8); % unexpanded flexible planning dam size [MCM]
        optParam.numPlan = x(9); % MCM
        optParam.planIncr = x(10);
        costParam.PercPlan = x(11); % initial upfront capital cost increase (0);
        costParam.PercPlanExp = x(12); % expansion cost of flexibly planned dam (0.5)
        costParam.discountrate = x(15);
        
        s_C = 1:3+optParam.numFlex+optParam.numPlan;
        M_C = length(s_C);
        
        storage = zeros(1, M_C);
        storage(1) = optParam.staticCap;
        storage(2) = optParam.smallFlexCap;
        storage(3) = optParam.smallPlanCap;
        storage(4:3+optParam.numFlex) = min(storage(2) + (1:optParam.numFlex)*optParam.flexIncr, 150);
        storage(4+optParam.numFlex:end) = min(storage(3) + (1:optParam.numPlan)*optParam.planIncr, 150);
        
        % Actions: Choose dam option in time period 1; 
        % expand dam in future time periods
        a_exp = 0:3+optParam.numFlex+optParam.numPlan;
        
        % Define infrastructure costs
        infra_cost = zeros(1,length(a_exp));
        infra_cost(2) = storage2damcost(storage(1),0); % cost of static dam
        for i = 1:optParam.numFlex % cost of flexible design dam
            [infra_cost(3), infra_cost(i+4)] = storage2damcost(storage(2), ...
                storage(i+3),costParam.PercFlex, costParam.PercFlexExp); % cost of flexible design exp to option X
        end
        for i = 1:optParam.numPlan % cost of flexible plan dam
            [infra_cost(4), infra_cost(i+4+optParam.numFlex)] = storage2damcost(storage(3), ...
                storage(i+3+optParam.numFlex),costParam.PercPlan, costParam.PercPlanExp); % cost of flexible planning exp to option X
        end
        
        if s == 1
            infra_cost_nonadaptive = struct();
            infra_cost_nonadaptive.static = [storage(1); infra_cost(2)/1E6];
            infra_cost_nonadaptive.flex = [[storage(2), storage(4:3+optParam.numFlex)]; [infra_cost(3), infra_cost(5:optParam.numFlex+4)]./1E6];
            infra_cost_nonadaptive.plan = [[storage(3), storage(4+optParam.numFlex:end)]; [infra_cost(4), infra_cost(5+optParam.numFlex:end)]./1E6];
            infra_cost_nonadaptive_lookup = infra_cost;
        else
            infra_cost_adaptive = struct();
            infra_cost_adaptive.static = [storage(1); infra_cost(2)/1E6];
            infra_cost_adaptive.flex = [[storage(2), storage(4:3+optParam.numFlex)]; [infra_cost(3), infra_cost(5:optParam.numFlex+4)]./1E6];
            infra_cost_adaptive.plan = [[storage(3), storage(4+optParam.numFlex:end)]; [infra_cost(4), infra_cost(5+optParam.numFlex:end)]./1E6];
            infra_cost_adaptive_lookup = infra_cost;
        end
    end

% =========================================================================

    % create bar plot of expected dam and shortage costs:
    f=figure('units','centimeters','position',[0,0,19,9]);
    t = tiledlayout(3,1);
    t.TileSpacing = 'tight';
    
    P_ranges = {[66:1:76];[77:1:86];[87:97]}; % reference dry, moderate, and wet climates
    
    s_C_adapt = [bestAct_adapt(2) bestAct_adapt(3), bestAct_adapt(8) (bestAct_adapt(3)+[1:bestAct_adapt(4)]*bestAct_adapt(5)),...
        (bestAct_adapt(8)+[1:bestAct_adapt(9)]*bestAct_adapt(10))];
    
    s_C_nonadapt = [bestAct_nonadapt(2) bestAct_nonadapt(3), bestAct_nonadapt(8) (bestAct_nonadapt(3)+[1:bestAct_nonadapt(4)]*bestAct_nonadapt(5)),...
        (bestAct_nonadapt(8)+[1:bestAct_nonadapt(9)]*bestAct_adapt(10))];
    
    % forward simulations of shortage and infrastructure costs
    totalCost_adapt = squeeze(totalCostTime_adapt(:,:,1:3))/1E6;% 1 is flex, 2 is static
    totalCost_nonadapt = squeeze(totalCostTime_nonadapt(:,:,1:3))/1E6;% 1 is flex, 2 is static
    damCost_adapt = squeeze(damCostTime_adapt(:,:,1:3))/1E6;% 1 is flex, 2 is static
    damCost_nonadapt = squeeze(damCostTime_nonadapt(:,:,1:3))/1E6;% 1 is flex, 2 is static
    
    for p = 1:length(P_ranges)
        % re-initialize mean cost arrays:
        totalCost_P = NaN(6,1);
        damCost_P = NaN(6,1);
        expCost_P = NaN(6,1);
        shortageCost_P = NaN(6,1);
        
        % calculate the average shortage, dam, and expansion costs for each
        % climate state:
        ind_P = ismember(P_state_adapt(:,5),P_ranges{p});
        
        % find climate subset from data for adaptive and non-adaptive
        % operations
        totalCost_P([2,1,3]) = squeeze(mean(sum(totalCost_nonadapt(ind_P,:,:),2)));
        damCost_P([2,1,3]) = squeeze(mean(damCost_nonadapt(ind_P,1,:)));
        expCost_P([2,1,3]) = squeeze(mean(sum(damCost_nonadapt(ind_P,2:5,:),2)));
        shortageCost_P([2,1,3]) = squeeze(mean(sum(totalCost_nonadapt(ind_P,2:5,:)- damCost_nonadapt(ind_P,2:5,:),2)));    
        
        totalCost_P([5,4,6]) = squeeze(mean(sum(totalCost_adapt(ind_P,:,:),2)));
        damCost_P([5,4,6]) = squeeze(mean(damCost_adapt(ind_P,1,:)));
        expCost_P([5,4,6]) = squeeze(mean(sum(damCost_adapt(ind_P,2:5,:),2)));
        shortageCost_P([5,4,6]) = squeeze(mean(sum(totalCost_adapt(ind_P,2:5,:)- damCost_adapt(ind_P,2:5,:),2)));    

        costs = [damCost_P, expCost_P, shortageCost_P];
        
        % bar plot of average costs from simulations:
        nexttile
        c = [[230,230,230];[110, 110, 110];[177,89,89]; [114,114,114];[90,90,90];[44,44,44]]/255;
        b = bar([1,2,3,4.2,5.2,6.2],costs,'stacked','FaceColor',"flat");
        if p == 1
        title({strcat("Dry Final Climates")},'FontSize',13)
        set(gca, 'XTick',[],'FontSize',0.1)
        h=gca; 
        h.XAxis.TickLength = [0 0];
        set(gca,'XTickLabel',[],'FontSize',0.1)
        leg = legend('Initial Dam Cost','Dam Expansion Cost','Shortage Cost','FontSize',6);
        set(leg,'Box','off')
        else
            title({strcat("Wet Final Climates")},'FontSize',13)
            if p==2
                title({strcat("Moderate Final Climates")},'FontSize',13)
                ylabel('Mean Simulated Cost (M$)','FontWeight','bold')
                set(gca, 'XTick',[],'FontSize',0.1)
                h=gca; 
                h.XAxis.TickLength = [0 0];
                set(gca,'XTickLabel',[],'FontSize',0.1)
            end        
        end

        for k = 1:size(b,2)
            b(k).CData = c(k,:);
        end
        
        ax = gca;
        ax.LineWidth = 0.5;

        xlim([0.5,6.7])
        yl = ylim;
        ylim([0, yl(2)+0.2*yl(2)])

        if p == 3

            set(gca, 'XTick', [1,2,3,4.2,5.2,6.2])
    
            set(gca,'XTickLabel',{strcat('Static Design'),...
                strcat('Flexible Design'),...
                strcat('Flexible Planning'),...
                strcat('Static Design'),...
                strcat('Flexible Design'),...
                strcat('Flexible Planning')})
            
            text(1.5, -35, "Static Operations",'FontSize',15,'FontWeight','bold')
            text(4.45, -35, "Flexible Operations",'FontSize',15,'FontWeight','bold')
                  
        end
    end
    
end

% save figure as .pdf and .eps and update fontsize
set(findall(gcf,'-property','FontSize'),'FontSize',8)
set(leg,'FontSize',6.5,'Position', leg.Position + [0.04 0.05 0 0])
cd(strcat(dir, 'Plots'))
exportgraphics(gcf, 'Figures/mainbody/fig6_barplot.pdf', 'ContentType', 'vector'); 
saveas(gcf, 'Figures/mainbody/fig6_barplot.eps')
saveas(gcf, 'Figures/mainbody/fig6_barplot.jpg')


%% Figure 7: Discount rate sensitivity plot

% Description: (a) Comparison of mean simulated discounted total costs
% across infrastructure alternatives under 0%, 3%, and 6% discounting for 
% simulations that transition to representative dry final climate states 
% (66-76 mm/mo). (b) Line plot of mean total costs vs. time for the subset 
% of 10,000 forward simulations that transition to representative dry 
% climate states (66-76 mm/mo) under 0%, 3%, and 6% discounting scenarios 
% for the flexible planning with flexible operations infrastructure 
% alternative. 

% =========================================================================
% (1) LOAD THE DATA FILES TO PLOT

discounts = [0, 3, 6];
cprimes = 1.25e-6;

% specify file path for data
cd('/Users/keaniw/Documents/Research/Kenya Project/Project Code/Flexible_Ops_Plan_Design/')

% initialize figure graphics
f=figure('units','centimeters','position',[0,0,21,11]);
t = tiledlayout(3,7);
t.TileSpacing = 'tight';
t.Padding = 'loose';

% ===========PANEL A: STACKED BAR FOR MEAN DISCOUNTED COSTS================

for d=1:length(discounts)
    disc = string(discounts(d));
    for c = 1:length(cprimes)
        cp = cprimes(c);
        c_prime = regexprep(strrep(string(cp), '.', ''), {'-0'}, {''});
        
        % specify folder
        folder = 'Results/Results_SDP_expansion';
        
        % load adaptive operations files:
        load(strcat(folder,'/BestFlex_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestStatic_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestPlan_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestFlexStaticPlan_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))

        bestAct_adapt = bestAct;
        V_adapt = V;
        Vs_adapt = Vs;
        Vd_adapt = Vd;
        Vs_static_adapt = allVs_static;
        Vs_flex_adapt = allVs_flex;
        Vs_plan_adapt = allVs_plan;
        Vd_static_adapt = allVd_static;
        Vd_flex_adapt = allVd_flex;
        Vd_plan_adapt = allVd_plan;
        X_adapt = X;
        C_adapt = C_state;
        action_adapt = action;
        totalCostTime_adapt = totalCostTime;
        damCostTime_adapt = damCostTime;
        P_state_adapt = P_state;
        bestVal_flex_adapt = bestVal_flex;
        bestVal_static_adapt = bestVal_static;
        bestVal_plan_adapt = bestVal_plan;
        if bestAct(3) + bestAct(4)*bestAct(5) > 150
            bestAct_adapt(4) = (150 - bestAct(2))/bestAct(5);
        end
        if bestAct(8) + bestAct(9)*bestAct(10) > 150
            bestAct_adapt(9) = (150 - bestAct(8))/bestAct(10);
        end
        
        % load non-adaptive operations files:
        load(strcat(folder,'/BestFlex_nonadaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestStatic_nonadaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestPlan_nonadaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestFlexStaticPlan_nonadaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))

        bestAct_nonadapt = bestAct;
        V_nonadapt = V;
        Vs_nonadapt = Vs;
        Vd_nonadapt = Vd;
        Vs_static_nonadapt = allVs_static;
        Vs_flex_nonadapt = allVs_flex;
        Vs_plan_nonadapt = allVs_plan;
        Vd_static_nonadapt = allVd_static;
        Vd_flex_nonadapt = allVd_flex;
        Vd_plan_nonadapt = allVd_plan;
        X_nonadapt = X;
        C_nonadapt = C_state;
        action_nonadapt = action;
        totalCostTime_nonadapt = totalCostTime;
        damCostTime_nonadapt = damCostTime;
        P_state_nonadapt = P_state;
        bestVal_flex_nonadapt = bestVal_flex;
        bestVal_static_nonadapt = bestVal_static;
        bestVal_plan_nonadapt = bestVal_plan;
        if bestAct(3) + bestAct(4)*bestAct(5) > 150
            bestAct_nonadapt(4) = (150 - bestAct(2))/bestAct(5);
        end
        if bestAct(8) + bestAct(9)*bestAct(10) > 150
            bestAct_nonadapt(9) = (150 - bestAct(8))/bestAct(10);
        end
        
    end
    
    % set different decade label parameters for use in plotting
    decade = {'2001-2020', '2021-2040', '2041-2060', '2061-2080', '2081-2100'};
    decade_short = {'2001-20', '2021-40', '2041-60', '2061-80', '2081-00'};
    decadeline = {'2001-\newline2020', '2021-\newline2040', '2041-\newline2060', '2061-\newline2080', '2081-\newline2100'};

    for s = 1:2 % for non-adaptive and adaptive operations
        if s == 1
            x = bestAct_nonadapt;
        else
            x = bestAct_adapt;
        end
        
        % load optimal infrastructure information
        optParam.staticCap = x(2); % static dam size [MCM]
        optParam.smallFlexCap = x(3); % unexpanded flexible design dam size [MCM]
        optParam.numFlex = x(4);  % number of possible expansion capacities [#]
        optParam.flexIncr = x(5); % increment of flexible expansion capacities [MCM]
        costParam.PercFlex = x(6); % Initial upfront capital cost increase (0.075)
        costParam.PercFlexExp = x(7); % Expansion cost of flexible dam  (0.15)
        optParam.smallPlanCap = x(8); % unexpanded flexible planning dam size [MCM]
        optParam.numPlan = x(9); % MCM
        optParam.planIncr = x(10);
        costParam.PercPlan = x(11); % initial upfront capital cost increase (0);
        costParam.PercPlanExp = x(12); % expansion cost of flexibly planned dam (0.5)
        costParam.discountrate = x(15);
        
        s_C = 1:3+optParam.numFlex+optParam.numPlan;
        M_C = length(s_C);
        
        storage = zeros(1, M_C);
        storage(1) = optParam.staticCap;
        storage(2) = optParam.smallFlexCap;
        storage(3) = optParam.smallPlanCap;
        storage(4:3+optParam.numFlex) = min(storage(2) + (1:optParam.numFlex)*optParam.flexIncr, 150);
        storage(4+optParam.numFlex:end) = min(storage(3) + (1:optParam.numPlan)*optParam.planIncr, 150);
        
        % Actions: Choose dam option in time period 1; expand dam in future time periods
        a_exp = 0:3+optParam.numFlex+optParam.numPlan;
        
        % Define infrastructure costs
        infra_cost = zeros(1,length(a_exp));
        infra_cost(2) = storage2damcost(storage(1),0); % cost of static dam
        for i = 1:optParam.numFlex % cost of flexible design dam
            [infra_cost(3), infra_cost(i+4)] = storage2damcost(storage(2), ...
                storage(i+3),costParam.PercFlex, costParam.PercFlexExp); % cost of flexible design exp to option X
        end
        for i = 1:optParam.numPlan % cost of flexible plan dam
            [infra_cost(4), infra_cost(i+4+optParam.numFlex)] = storage2damcost(storage(3), ...
                storage(i+3+optParam.numFlex),costParam.PercPlan, costParam.PercPlanExp); % cost of flexible planning exp to option X
        end
        
        if s == 1
            infra_cost_nonadaptive = struct();
            infra_cost_nonadaptive.static = [storage(1); infra_cost(2)/1E6];
            infra_cost_nonadaptive.flex = [[storage(2), storage(4:3+optParam.numFlex)]; [infra_cost(3), infra_cost(5:optParam.numFlex+4)]./1E6];
            infra_cost_nonadaptive.plan = [[storage(3), storage(4+optParam.numFlex:end)]; [infra_cost(4), infra_cost(5+optParam.numFlex:end)]./1E6];
            infra_cost_nonadaptive_lookup = infra_cost/1E6;
        else
            infra_cost_adaptive = struct();
            infra_cost_adaptive.static = [storage(1); infra_cost(2)/1E6];
            infra_cost_adaptive.flex = [[storage(2), storage(4:3+optParam.numFlex)]; [infra_cost(3), infra_cost(5:optParam.numFlex+4)]./1E6];
            infra_cost_adaptive.plan = [[storage(3), storage(4+optParam.numFlex:end)]; [infra_cost(4), infra_cost(5+optParam.numFlex:end)]./1E6];
            infra_cost_adaptive_lookup = infra_cost/1E6;
        end
    end
    
    facecolors = [[153,204,204]; [204,255,255]; [153, 153, 204]; [204, 204, 255];...
        [255, 102, 102]; [255, 153, 153]]/255;
    
    P_regret = [{66:76}]; % dry, moderate, wet
    
    s_C_adapt = [bestAct_adapt(2) bestAct_adapt(3), bestAct_adapt(8) (bestAct_adapt(3)+[1:bestAct_adapt(4)]*bestAct_adapt(5)),...
        (bestAct_adapt(8)+[1:bestAct_adapt(9)]*bestAct_adapt(10))];
    
    s_C_nonadapt = [bestAct_nonadapt(2) bestAct_nonadapt(3), bestAct_nonadapt(8) (bestAct_nonadapt(3)+[1:bestAct_nonadapt(4)]*bestAct_nonadapt(5)),...
        (bestAct_nonadapt(8)+[1:bestAct_nonadapt(9)]*bestAct_adapt(10))];
    
    P_ranges = {[66:1:76]}; % reference dry, moderate, and wet climates

    
    % forward simulations of shortage and infrastructure costs
    totalCost_adapt = squeeze(totalCostTime_adapt(:,:,1:3))/1E6;% 1 is flex, 2 is static
    totalCost_nonadapt = squeeze(totalCostTime_nonadapt(:,:,1:3))/1E6;% 1 is flex, 2 is static
    damCost_adapt = squeeze(damCostTime_adapt(:,:,1:3))/1E6;% 1 is flex, 2 is static
    damCost_nonadapt = squeeze(damCostTime_nonadapt(:,:,1:3))/1E6;% 1 is flex, 2 is static
    
    for p = 1:length(P_ranges)
        
        % re-initialize mean cost arrays:
        totalCost_P = NaN(6,1);
        damCost_P = NaN(6,1);
        expCost_P = NaN(6,1);
        shortageCost_P = NaN(6,1);
        
        % calculate the average shortage, dam, and expansion costs for each
        % climate state:
        ind_P = ismember(P_state_adapt(:,5),P_ranges{p});
        
        % find climate subset from data for adaptive and non-adaptive
        % operations
        totalCost_P([2,1,3]) = squeeze(mean(sum(totalCost_nonadapt(ind_P,:,:),2)));
        damCost_P([2,1,3]) = squeeze(mean(damCost_nonadapt(ind_P,1,:)));
        expCost_P([2,1,3]) = squeeze(mean(sum(damCost_nonadapt(ind_P,2:5,:),2)));
        shortageCost_P([2,1,3]) = squeeze(mean(sum(totalCost_nonadapt(ind_P,2:5,:)- damCost_nonadapt(ind_P,2:5,:),2)));    
        
        totalCost_P([5,4,6]) = squeeze(mean(sum(totalCost_adapt(ind_P,:,:),2)));
        damCost_P([5,4,6]) = squeeze(mean(damCost_adapt(ind_P,1,:)));
        expCost_P([5,4,6]) = squeeze(mean(sum(damCost_adapt(ind_P,2:5,:),2)));
        shortageCost_P([5,4,6]) = squeeze(mean(sum(totalCost_adapt(ind_P,2:5,:)- damCost_adapt(ind_P,2:5,:),2)));    

        costs = [damCost_P, expCost_P, shortageCost_P];
        
        % bar plot of average costs from simulations:
        nexttile([1,4])

        % set background color annotation
        patch([5.7 6.7 6.7 5.7], [400 400 0 0], [0.92 0.92 0.92],'EdgeColor', ...
            'none')
        hold on

        c = [[230,230,230];[110, 110, 110];[177,89,89]; [114,114,114];[90,90,90];[44,44,44]]/255;
        b = bar([1,2,3,4.2,5.2,6.2],costs,'stacked','FaceColor',"flat");
        if d == 3
            ax = gca;
            ax.XTick = [1,2,3,4.2,5.2,6.2];
            ax.XTickLabel = '';
            myLabels = { 'Static', 'Flexible', 'Flexible', 'Static', 'Flexible', 'Flexible'; 
                'Design', 'Design', 'Planning', 'Design', 'Design', 'Planning'; 
                ' ', ' ', ' ',' ', ' ', ' '; ' ', ' ', ' ',' ', ' ', ' '};
            for i = 1:length(myLabels)
                text(ax.XTick(i), ax.YLim(1)-5, sprintf('%s\n%s\n%s\n%s\n', myLabels{:,i}), ...
                    'horizontalalignment', 'center', 'verticalalignment', 'top');    
            end
        
            text(1.2, -150, "Static Operations",'FontSize',8,'FontWeight','bold')
            text(4.1, -150, "Flexible Operations",'FontSize',8,'FontWeight','bold')
                  
            l = legend(b, 'Initial Dam Cost','Dam Expansion Cost','Shortage Cost','FontSize',8);
            set(l,'Box','off','Location','northwest')
            l.ItemTokenSize = [20,18];

          
            
        else
            set(gca, 'XTick',[],'FontSize',8)  
            if d==2
                ylabel('Discounted Cost (M$)','FontWeight','bold','FontSize',8)
            end
            if d==1
                % panel label
                text(-0.2, 420, "(a)",'FontSize',12,'FontWeight','bold')
            end
        end
        
        for k = 1:size(b,2)
            b(k).CData = c(k,:);
        end
        
        set(gca,'Box','on');
        set(gca, 'Layer', 'top')
        ax = gca;
         ax.LineWidth = 0.5;

        xlim([0.5,6.7])
        
        yl = ylim;
        ylim([0, 400])

       
    end
    title(strcat("Discount Rate: ", disc,'%'),'FontWeight','bold','FontSize',11)
    set(findall(gcf,'-property','FontSize'),'FontSize',8)

end
    

% ===============PANEL B: LINE PLOT OF NON-DISCOUNTED COSTS================

for d=1:length(discounts)
    disc = string(discounts(d));
    for c = 1:length(cprimes)
        cp = cprimes(c);
        c_prime = regexprep(strrep(string(cp), '.', ''), {'-0'}, {''});
        
        % re-specify folder
        folder = 'Results/Results_SDP_expansion';
        
        % load adaptive operations files:
        load(strcat(folder,'/BestFlex_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestStatic_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestPlan_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestFlexStaticPlan_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))

        bestAct_adapt = bestAct;
        V_adapt = V;
        Vs_adapt = Vs;
        Vd_adapt = Vd;
        Vs_static_adapt = allVs_static;
        Vs_flex_adapt = allVs_flex;
        Vs_plan_adapt = allVs_plan;
        Vd_static_adapt = allVd_static;
        Vd_flex_adapt = allVd_flex;
        Vd_plan_adapt = allVd_plan;
        X_adapt = X;
        C_adapt = C_state;
        action_adapt = action;
        totalCostTime_adapt = totalCostTime;
        damCostTime_adapt = damCostTime;
        P_state_adapt = P_state;
        bestVal_flex_adapt = bestVal_flex;
        bestVal_static_adapt = bestVal_static;
        bestVal_plan_adapt = bestVal_plan;
        if bestAct(3) + bestAct(4)*bestAct(5) > 150
            bestAct_adapt(4) = (150 - bestAct(2))/bestAct(5);
        end
        if bestAct(8) + bestAct(9)*bestAct(10) > 150
            bestAct_adapt(9) = (150 - bestAct(8))/bestAct(10);
        end
        
        % load non-adaptive operations files:
        load(strcat(folder,'/BestFlex_nonadaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestStatic_nonadaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestPlan_nonadaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestFlexStaticPlan_nonadaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))

        bestAct_nonadapt = bestAct;
        V_nonadapt = V;
        Vs_nonadapt = Vs;
        Vd_nonadapt = Vd;
        Vs_static_nonadapt = allVs_static;
        Vs_flex_nonadapt = allVs_flex;
        Vs_plan_nonadapt = allVs_plan;
        Vd_static_nonadapt = allVd_static;
        Vd_flex_nonadapt = allVd_flex;
        Vd_plan_nonadapt = allVd_plan;
        X_nonadapt = X;
        C_nonadapt = C_state;
        action_nonadapt = action;
        totalCostTime_nonadapt = totalCostTime;
        damCostTime_nonadapt = damCostTime;
        P_state_nonadapt = P_state;
        bestVal_flex_nonadapt = bestVal_flex;
        bestVal_static_nonadapt = bestVal_static;
        bestVal_plan_nonadapt = bestVal_plan;
        if bestAct(3) + bestAct(4)*bestAct(5) > 150
            bestAct_nonadapt(4) = (150 - bestAct(2))/bestAct(5);
        end
        if bestAct(8) + bestAct(9)*bestAct(10) > 150
            bestAct_nonadapt(9) = (150 - bestAct(8))/bestAct(10);
        end
        
    end
    
    % set different decade label parameters for use in plotting
    decade = {'2001-2020', '2021-2040', '2041-2060', '2061-2080', '2081-2100'};
    decade_short = {'2001-20', '2021-40', '2041-60', '2061-80', '2081-00'};
    decadeline = {'2001-\newline2020', '2021-\newline2040', '2041-\newline2060', '2061-\newline2080', '2081-\newline2100'};

        for s = 1:2 % for non-adaptive and adaptive operations
        if s == 1
            x = bestAct_nonadapt;
        else
            x = bestAct_adapt;
        end
        
        % load optimal infrastructure information
        optParam.staticCap = x(2); % static dam size [MCM]
        optParam.smallFlexCap = x(3); % unexpanded flexible design dam size [MCM]
        optParam.numFlex = x(4);  % number of possible expansion capacities [#]
        optParam.flexIncr = x(5); % increment of flexible expansion capacities [MCM]
        costParam.PercFlex = x(6); % Initial upfront capital cost increase (0.075)
        costParam.PercFlexExp = x(7); % Expansion cost of flexible dam  (0.15)
        optParam.smallPlanCap = x(8); % unexpanded flexible planning dam size [MCM]
        optParam.numPlan = x(9); % MCM
        optParam.planIncr = x(10);
        costParam.PercPlan = x(11); % initial upfront capital cost increase (0);
        costParam.PercPlanExp = x(12); % expansion cost of flexibly planned dam (0.5)
        costParam.discountrate = x(15);
        
        s_C = 1:3+optParam.numFlex+optParam.numPlan;
        M_C = length(s_C);
        
        storage = zeros(1, M_C);
        storage(1) = optParam.staticCap;
        storage(2) = optParam.smallFlexCap;
        storage(3) = optParam.smallPlanCap;
        storage(4:3+optParam.numFlex) = min(storage(2) + (1:optParam.numFlex)*optParam.flexIncr, 150);
        storage(4+optParam.numFlex:end) = min(storage(3) + (1:optParam.numPlan)*optParam.planIncr, 150);
        
        % Actions: Choose dam option in time period 1; expand dam in future time periods
        a_exp = 0:3+optParam.numFlex+optParam.numPlan;
        
        % Define infrastructure costs
        infra_cost = zeros(1,length(a_exp));
        infra_cost(2) = storage2damcost(storage(1),0); % cost of static dam
        for i = 1:optParam.numFlex % cost of flexible design dam
            [infra_cost(3), infra_cost(i+4)] = storage2damcost(storage(2), ...
                storage(i+3),costParam.PercFlex, costParam.PercFlexExp); % cost of flexible design exp to option X
        end
        for i = 1:optParam.numPlan % cost of flexible plan dam
            [infra_cost(4), infra_cost(i+4+optParam.numFlex)] = storage2damcost(storage(3), ...
                storage(i+3+optParam.numFlex),costParam.PercPlan, costParam.PercPlanExp); % cost of flexible planning exp to option X
        end
        
        if s == 1
            infra_cost_nonadaptive = struct();
            infra_cost_nonadaptive.static = [storage(1); infra_cost(2)/1E6];
            infra_cost_nonadaptive.flex = [[storage(2), storage(4:3+optParam.numFlex)]; [infra_cost(3), infra_cost(5:optParam.numFlex+4)]./1E6];
            infra_cost_nonadaptive.plan = [[storage(3), storage(4+optParam.numFlex:end)]; [infra_cost(4), infra_cost(5+optParam.numFlex:end)]./1E6];
            infra_cost_nonadaptive_lookup = infra_cost/1E6;
        else
            infra_cost_adaptive = struct();
            infra_cost_adaptive.static = [storage(1); infra_cost(2)/1E6];
            infra_cost_adaptive.flex = [[storage(2), storage(4:3+optParam.numFlex)]; [infra_cost(3), infra_cost(5:optParam.numFlex+4)]./1E6];
            infra_cost_adaptive.plan = [[storage(3), storage(4+optParam.numFlex:end)]; [infra_cost(4), infra_cost(5+optParam.numFlex:end)]./1E6];
            infra_cost_adaptive_lookup = infra_cost/1E6;
        end
        end
    
    % Non-discounted Boxplot of costs over time for dry climates
    
    facecolors = [[153,204,204]; [204,255,255]; [153, 153, 204]; [204, 204, 255];...
        [255, 102, 102]; [255, 153, 153]]/255;
    
    P_regret = [{66:76}]; % dry, moderate, wet
    
    s_C_adapt = [bestAct_adapt(2) bestAct_adapt(3), bestAct_adapt(8) (bestAct_adapt(3)+[1:bestAct_adapt(4)]*bestAct_adapt(5)),...
        (bestAct_adapt(8)+[1:bestAct_adapt(9)]*bestAct_adapt(10))];
    
    s_C_nonadapt = [bestAct_nonadapt(2) bestAct_nonadapt(3), bestAct_nonadapt(8) (bestAct_nonadapt(3)+[1:bestAct_nonadapt(4)]*bestAct_nonadapt(5)),...
        (bestAct_nonadapt(8)+[1:bestAct_nonadapt(9)]*bestAct_adapt(10))];
    
    % forward simulations of shortage and infrastructure costs
    totalCost_adapt = squeeze(totalCostTime_adapt(:,:,1:3))/1E6;% 1 is flex, 2 is static
    totalCost_nonadapt = squeeze(totalCostTime_nonadapt(:,:,1:3))/1E6;% 1 is flex, 2 is static
    damCost_adapt = squeeze(damCostTime_adapt(:,:,1:3))/1E6;% 1 is flex, 2 is static
    damCost_nonadapt = squeeze(damCostTime_nonadapt(:,:,1:3))/1E6;% 1 is flex, 2 is static
    
    for p = 1:length(P_regret) % for each transition to drier or wetter climate
        
        % Index realizations where final P state is 72, 79, or 87 mm/mo
        [ind_P_adapt,~] = find(P_state_adapt(:,end) == P_regret{1});
        [ind_P_nonadapt,~] = find(P_state_nonadapt(:,end) == P_regret{1});
        
        % P state time series: select associated time series
        Pnow_adapt = P_state_adapt(ind_P_adapt,:);
        Pnow_nonadapt = P_state_adapt(ind_P_nonadapt,:);
        
        % Actions: select associated simulations
        ActionPnow_adapt = action_adapt(ind_P_adapt,:,:);
        ActionPnow_nonadapt = action_nonadapt(ind_P_nonadapt,:,:);
        
        % map actions to dam capacity of time (fill actions with value 0)
        % == FLEXIBLE DESIGN ==
        actionPnowFlex_adapt = ActionPnow_adapt(:,:,1);
        for r = 1:size(actionPnowFlex_adapt,1) % each forward simulation
            for ia = 2:size(actionPnowFlex_adapt,2) % each subsequent time period
                if actionPnowFlex_adapt(r,ia) == 0
                    actionPnowFlex_adapt(r,ia) = actionPnowFlex_adapt(r,ia-1);
                end
            end
        end
        
        actionPnowFlex_nonadapt = ActionPnow_nonadapt(:,:,1);
        for r = 1:size(actionPnowFlex_nonadapt,1)
            for ia = 2:size(actionPnowFlex_nonadapt,2)
                if actionPnowFlex_nonadapt(r,ia) == 0
                    actionPnowFlex_nonadapt(r,ia) = actionPnowFlex_nonadapt(r,ia-1);
                end
            end
        end
        
        % == FLEXIBLE PLANNING ==
        actionPnowPlan_adapt = ActionPnow_adapt(:,:,3);
        for r = 1:size(actionPnowPlan_adapt,1)
            for ia = 2:size(actionPnowPlan_adapt,2)
                if actionPnowPlan_adapt(r,ia) == 0
                    actionPnowPlan_adapt(r,ia) = actionPnowPlan_adapt(r,ia-1);
                end
            end
        end
        
        actionPnowPlan_nonadapt = ActionPnow_nonadapt(:,:,3);
        for r = 1:size(actionPnowPlan_nonadapt,1)
            for ia = 2:size(actionPnowPlan_nonadapt,2)
                if actionPnowPlan_nonadapt(r,ia) == 0
                    actionPnowPlan_nonadapt(r,ia) = actionPnowPlan_nonadapt(r,ia-1);
                end
            end
        end
        
        % == STATIC DAM ==
        actionPnowStatic_adapt = ActionPnow_adapt(:,:,2);
        for r = 1:size(actionPnowStatic_adapt,1) % each forward simulation
            for ia = 2:size(actionPnowStatic_adapt,2) % each subsequent time period
                if actionPnowStatic_adapt(r,ia) == 0
                    actionPnowStatic_adapt(r,ia) = actionPnowStatic_adapt(r,ia-1);
                end
            end
        end
        
        actionPnowStatic_nonadapt = ActionPnow_nonadapt(:,:,2);
        for r = 1:size(actionPnowStatic_nonadapt,1)
            for ia = 2:size(actionPnowStatic_nonadapt,2)
                if actionPnowStatic_nonadapt(r,ia) == 0
                    actionPnowStatic_nonadapt(r,ia) = actionPnowStatic_nonadapt(r,ia-1);
                end
            end
        end
        
        
        if disc == "0" % use discounted costs
            % Dam Cost Time: select associated simulations for given P state
            damCostPnowFlex_adapt = damCost_adapt(ind_P_adapt,:,1);
            damCostPnowFlex_nonadapt = damCost_nonadapt(ind_P_adapt,:,1);
            damCostPnowStatic_adapt = damCost_adapt(ind_P_adapt,:,2);
            damCostPnowStatic_nonadapt = damCost_nonadapt(ind_P_adapt,:,2);
            damCostPnowPlan_adapt = damCost_adapt(ind_P_adapt,:,3);
            damCostPnowPlan_nonadapt = damCost_nonadapt(ind_P_adapt,:,3);
            
            % Shortage Cost Time: select associated simulations
            shortageCostPnowFlex_adapt = totalCost_adapt(ind_P_adapt,:,1)-damCost_adapt(ind_P_adapt,:,1);
            shortageCostPnowFlex_nonadapt = totalCost_nonadapt(ind_P_adapt,:,1)-damCost_nonadapt(ind_P_adapt,:,1);
            shortageCostPnowStatic_adapt = totalCost_adapt(ind_P_adapt,:,2)-damCost_adapt(ind_P_adapt,:,2);
            shortageCostPnowStatic_nonadapt = totalCost_nonadapt(ind_P_adapt,:,2)-damCost_nonadapt(ind_P_adapt,:,2);
            shortageCostPnowPlan_adapt = totalCost_adapt(ind_P_adapt,:,3)-damCost_adapt(ind_P_adapt,:,3);
            shortageCostPnowPlan_nonadapt = totalCost_nonadapt(ind_P_adapt,:,3)-damCost_nonadapt(ind_P_adapt,:,3);
            
        else % non-discounted:use shortage cost files
            for i=1:5  
                
                % dam capacity at that time
                staticCap_adapt = s_C_adapt(actionPnowStatic_adapt(:,i));
                staticCap_nonadapt = s_C_nonadapt(actionPnowStatic_nonadapt(:,i));
                flexCap_adapt = s_C_adapt(actionPnowFlex_adapt(:,i));
                flexCap_nonadapt = s_C_nonadapt(actionPnowFlex_nonadapt(:,i));
                planCap_adapt = s_C_adapt(actionPnowPlan_adapt(:,i));
                planCap_nonadapt = s_C_nonadapt(actionPnowPlan_nonadapt(:,i)); % flex plan capacity at N=i
                
                % dam cost based on capacity at that time (infra_cost tables)
                
                % damCostPnowStatic_adapt = interp1(infra_cost_adaptive.static(1,:), infra_cost_adaptive.static(2,:), staticCap_adapt);
                damCostPnowStatic_adapt(:,i) = infra_cost_adaptive_lookup(ActionPnow_adapt(:,i,2)+1);
                damCostPnowStatic_nonadapt(:,i) = infra_cost_nonadaptive_lookup(ActionPnow_nonadapt(:,i,2)+1);
                damCostPnowFlex_adapt(:,i) = infra_cost_adaptive_lookup(ActionPnow_adapt(:,i,1)+1);
                damCostPnowFlex_nonadapt(:,i) = infra_cost_nonadaptive_lookup(ActionPnow_nonadapt(:,i,1)+1);
                damCostPnowPlan_adapt(:,i) = infra_cost_adaptive_lookup(ActionPnow_adapt(:,i,3)+1);
                damCostPnowPlan_nonadapt(:,i) = infra_cost_nonadaptive_lookup(ActionPnow_nonadapt(:,i,3)+1);
                
                % corresponding shortage cost at that time and P state
                if i==1
                    shortageCostPnowStatic_nonadapt(:,i)=0;
                    shortageCostPnowStatic_adapt(:,i)=0;
                    shortageCostPnowFlex_nonadapt(:,i)=0;
                    shortageCostPnowFlex_adapt(:,i)=0;
                    shortageCostPnowPlan_nonadapt(:,i)=0;
                    shortageCostPnowPlan_adapt(:,i)=0;
                else
                    folder = 'Results/Results_SDP_reservoir_ops';

                    for j=1:length(staticCap_nonadapt)
                        Ps = [66:1:97]; % P states indexed for 18:49 (Keani)
                        P_now = Pnow_adapt(j,i);
                        
                        s_state_filename = strcat(string(folder),'/sdp_nonadaptive_shortage_cost_s',string(staticCap_nonadapt(j)),'.mat');
                        shortageCostDir = load(s_state_filename,'shortageCost');
                        shortageCost_s_state = shortageCostDir.shortageCost(:,18:49)/1E6; % 66 mm/month to 97 mm/month
                        shortageCostPnowStatic_nonadapt(j,i) = shortageCost_s_state(i, Ps == P_now)*cp;
                        s_state_filename = strcat(string(folder),'/sdp_adaptive_shortage_cost_s',string(staticCap_adapt(j)),'.mat');
                        shortageCostDir = load(s_state_filename,'shortageCost');
                        shortageCost_s_state = shortageCostDir.shortageCost(:,18:49)/1E6; % 66 mm/month to 97 mm/month
                        shortageCostPnowStatic_adapt(j,i) = shortageCost_s_state(i, Ps == P_now)*cp;
                        
                        s_state_filename = strcat(string(folder),'/sdp_nonadaptive_shortage_cost_s',string(flexCap_nonadapt(j)),'.mat');
                        shortageCostDir = load(s_state_filename,'shortageCost');
                        shortageCost_s_state = shortageCostDir.shortageCost(:,18:49)/1E6; % 66 mm/month to 97 mm/month
                        shortageCostPnowFlex_nonadapt(j,i) = shortageCost_s_state(i, Ps == P_now)*cp;
                        s_state_filename = strcat(string(folder),'/sdp_adaptive_shortage_cost_s',string(flexCap_adapt(j)),'.mat');
                        shortageCostDir = load(s_state_filename,'shortageCost');
                        shortageCost_s_state = shortageCostDir.shortageCost(:,18:49)/1E6; % 66 mm/month to 97 mm/month
                        shortageCostPnowFlex_adapt(j,i) = shortageCost_s_state(i, Ps == P_now)*cp;
                        
                        s_state_filename = strcat(string(folder),'/sdp_nonadaptive_shortage_cost_s',string(planCap_nonadapt(j)),'.mat');
                        shortageCostDir = load(s_state_filename,'shortageCost');
                        shortageCost_s_state = shortageCostDir.shortageCost(:,18:49)/1E6; % 66 mm/month to 97 mm/month
                        shortageCostPnowPlan_nonadapt(j,i) = shortageCost_s_state(i, Ps == P_now)*cp;
                        s_state_filename = strcat(string(folder),'/sdp_adaptive_shortage_cost_s',string(planCap_adapt(j)),'.mat');
                        shortageCostDir = load(s_state_filename,'shortageCost');
                        shortageCost_s_state = shortageCostDir.shortageCost(:,18:49)/1E6; % 66 mm/month to 97 mm/month
                        shortageCostPnowPlan_adapt(j,i) = shortageCost_s_state(i, Ps == P_now)*cp;
                    end
                end
            end
        end
        
        % =================================================================
        % line plot mean non-discounted cost over time (only for flexible planning):
        nondisc_DamCostTime = damCostPnowPlan_adapt + shortageCostPnowPlan_adapt;
        cols = [[132 169 140]; [102, 155, 188]; [53 79 82]]/255;
        
        if d==1
            %tmerged = nexttile([3,1]);
            tmerged = nexttile([3,3]);
        else
            axes(tmerged)
        end
        hold on

        % set background color annotation
        if d==1
            patch([0.75 5.25 5.25 0.75], [100 100 0 0], [0.92 0.92 0.92], ...
                'EdgeColor', ...
                'none')
            hold on
        end

        p1(d)=plot((1:5),mean(nondisc_DamCostTime), '-','Marker','.', ...
            'MarkerSize',20,'LineWidth',2,...
            'Color',cols(d,:),'MarkerEdgeColor',cols(d,:));
        
        if d==3
            xticks([1 2 3 4 5])
            xticklabels(decade_short)
            xlabel('Time Period','Fontweight','bold')
            ylabel('\newline Non-Discounted Cost (M$)','FontWeight','bold')
            
            l = legend([p1(3) p1(2) p1(1)],{'6%','3%','0%'},'Location','southeast');
            set(l,'Box','off')
            title(l,'Discount Rate')
        end

        
        if d==3
            xticks([1, 2, 3, 4, 5])
            xticklabels(decade_short)
            xlabel('Time Period','Fontweight','bold')

            set(gca,'Box','on');
            set(gca, 'Layer', 'top')
            
            ax.YAxis.Color = 'k';
            ax.LineWidth = 0.5;

            % panel label
            text(-0.03, 102, "(b)",'FontSize',12,'FontWeight','bold')
        else
            xticks([0, 1, 2, 3, 4])
            xticklabels([])
        end
    end
       
    xlim([0.75, 5.25])
    ylim([0,100])
    
end


% save figure as .pdf and .eps and update fontsize
set(findall(gcf,'-property','FontSize'),'FontSize',8)
cd(strcat(dir, 'Plots'))
exportgraphics(gcf, 'Figures/mainbody/fig7_disc.pdf', 'ContentType', 'vector'); 
saveas(gcf, 'Figures/mainbody/fig7_disc.eps')
saveas(gcf, 'Figures/mainbody/fig7_disc.jpg')


%% Figure 8: Comparison of Mean Shortage Cost by Climate (c'=1.25E-6 vs 1E-7)

% Description: (a) Initial capacity of optimal dam designs selected by the 
% infrastructure screening model across a range of values of c: under 0% 
% discounting. (b) Comparison of simulated mean shortage costs for a range 
% of representative final climate states and for scenarios c’base and 
% c’alt. Results presented in Figure 8 assume 0% discounting.

% =========================================================================
% (1) LOAD THE DATA FILES TO PLOT

f=figure('units','centimeters','position',[0,0,19,12]);

% specify file path for data
cd('/Users/keaniw/Documents/Research/Kenya Project/Project Code/Flexible_Ops_Plan_Design/')

% =========================================================================
% (2) PANEL A: SCATTER OPTIMAL INITIAL DAM CAPACITY BY C'

% Make scatter plots of values
folder = 'Results_cprime_sensitivity'; % folder to load
facecolors = [[153,204,204]; [153, 153, 204]; ...
        [255, 102, 102]; [255, 153, 153]]/255;

subplot(1,7,[1:3])

for des = 1:3 % (1) static design (2) flexible design (3) flexible planning
    for ops = 1:2 % (1) static ops (2) flexible ops
        load(strcat(folder,'/cprimes_des',string(des),'_ops',string(ops),'_disc0.mat'))
        if ops == 1 % flexible ops
            s(des,ops)=scatter(round_min_cprimes,60:10:140,50,'o','LineWidth',1.5,...
                'MarkerEdgeColor',facecolors(des,:));
            if des == 1 % for static ops, add point to 150 MCM
                hold on
            end
        else
            s(des,ops)=scatter(round_min_cprimes,60:10:140,50,'*','LineWidth',1.5,...
                'MarkerEdgeColor',facecolors(des,:));
            if des == 1 % for static ops, add point to 150 MCM
                hold on
            end
        end
        hold on
    end
end
l1 = xline(1.25E-6,'--','DisplayName',"c'_{alt}",'LineWidth',1.5);
hold on
l2 = xline(1.00E-7,'-.','DisplayName',"c'_{base}",'LineWidth',1.5);
labels=["c'_{alt}","c'_{base}"];
t = text([1.00E-7+0.1E-6, 1.25E-6+0.1E-6],[138 138],labels,"HorizontalAlignment","left","VerticalAlignment","bottom",'FontSize',12);
set(t,'Rotation',0);

ylim([50,145])
xlim([1E-8, 3E-6])
xlabel("c' ($/m^6)",'FontWeight','bold')
ylabel("Optimal Initial Dam Capacity (MCM)\newline ",'FontWeight','bold')
box on
ax=gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;

hold on
pt(1)=scatter(-10,-10,'*','LineWidth',1.5,...
                'MarkerEdgeColor',[0 0 0]/255);
hold on
pt(2)=scatter(-10,-10,'o','LineWidth',1.5,...
                'MarkerEdgeColor',[0 0 0]);
hold on
dtype(1) = bar(-10,-10);
dtype(1).FaceColor = facecolors(1,:);
hold on
dtype(2) = bar(-10,-10);
dtype(2).FaceColor = facecolors(2,:);
hold on
dtype(3) = bar(-10,-10);
dtype(3).FaceColor = facecolors(3,:);

leg1 = legend([pt(1), pt(2)],'Static','Flexible');
leg1.FontSize = 10;
set(leg1,'Box','off','Location','southeast')
title(leg1,'Operations')

ah1=axes('position',get(gca,'position'),'visible','off');
leg2 = legend(ah1,[dtype(1) dtype(2) dtype(3)],'Static Design','Flexible Design', 'Flexible Planning');
leg2.FontSize = 10;
leg2.LineWidth = 1.5;
set(leg2,'Box','off')
set(leg2,'Location','southeast')
set(leg2,'Position',leg2.Position + [0.03 0 0 0])
set(leg1,'Position',leg2.Position + [0 0.15 0 0])
leg2.ItemTokenSize = [20,18];
leg1.ItemTokenSize = [20,18];
title(leg2, 'Design')

% panel label
text(-0.17, 1.02, "(a)",'FontSize',12,'FontWeight','bold')

hold off

% =========================================================================
% (3) PANEL B: SCATTER MEAN SHORTAGE COST BY CLIMATE FOR C'=1.25E-6 OR 1E-7


% specify file naming convensions
discounts = 0;
cprimes = [1.25e-6 1E-7];
            
subplot('Position',subplot(1,7,[4:7]).Position + [0.06 0 0 0])

% respecify folder
folder = 'Results_SDP_expansion'; % folder to load

for d=1
    disc = string(discounts(d));
    for cs = 1:2
        cp = cprimes(cs);
        c_prime = regexprep(strrep(string(cp), '.', ''), {'-0'}, {''});
            
        % load adaptive operations files:
        load(strcat(folder,'/BestFlex_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestStatic_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestPlan_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestFlexStaticPlan_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        
        bestAct_adapt = bestAct;
        V_adapt = V;
        Vs_adapt = Vs;
        Vd_adapt = Vd;
        Vs_static_adapt = allVs_static;
        Vs_flex_adapt = allVs_flex;
        Vs_plan_adapt = allVs_plan;
        Vd_static_adapt = allVd_static;
        Vd_flex_adapt = allVd_flex;
        Vd_plan_adapt = allVd_plan;
        X_adapt = X;
        C_adapt = C_state;
        action_adapt = action;
        totalCostTime_adapt = totalCostTime;
        damCostTime_adapt = damCostTime;
        P_state_adapt = P_state;
        bestVal_flex_adapt = bestVal_flex;
        bestVal_static_adapt = bestVal_static;
        bestVal_plan_adapt = bestVal_plan;
        if bestAct(3) + bestAct(4)*bestAct(5) > 150
            bestAct_adapt(4) = (150 - bestAct(2))/bestAct(5);
        end
        if bestAct(8) + bestAct(9)*bestAct(10) > 150
            bestAct_adapt(9) = (150 - bestAct(8))/bestAct(10);
        end
        
        % load non-adaptive operations files:
        load(strcat(folder,'/BestFlex_nonadaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestStatic_nonadaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestPlan_nonadaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
        load(strcat(folder,'/BestFlexStaticPlan_nonadaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))

        bestAct_nonadapt = bestAct;
        V_nonadapt = V;
        Vs_nonadapt = Vs;
        Vd_nonadapt = Vd;
        Vs_static_nonadapt = allVs_static;
        Vs_flex_nonadapt = allVs_flex;
        Vs_plan_nonadapt = allVs_plan;
        Vd_static_nonadapt = allVd_static;
        Vd_flex_nonadapt = allVd_flex;
        Vd_plan_nonadapt = allVd_plan;
        X_nonadapt = X;
        C_nonadapt = C_state;
        action_nonadapt = action;
        totalCostTime_nonadapt = totalCostTime;
        damCostTime_nonadapt = damCostTime;
        P_state_nonadapt = P_state;
        bestVal_flex_nonadapt = bestVal_flex;
        bestVal_static_nonadapt = bestVal_static;
        bestVal_plan_nonadapt = bestVal_plan;
        if bestAct(3) + bestAct(4)*bestAct(5) > 150
            bestAct_nonadapt(4) = (150 - bestAct(2))/bestAct(5);
        end
        if bestAct(8) + bestAct(9)*bestAct(10) > 150
            bestAct_nonadapt(9) = (150 - bestAct(8))/bestAct(10);
        end
        
        P_ranges = {[66:76],[77:86],[87:97]}; % reference dry, moderate, and wet climates
        
        s_C_adapt = [bestAct_adapt(2) bestAct_adapt(3), bestAct_adapt(8) (bestAct_adapt(3)+[1:bestAct_adapt(4)]*bestAct_adapt(5)),...
            (bestAct_adapt(8)+[1:bestAct_adapt(9)]*bestAct_adapt(10))];
        
        s_C_nonadapt = [bestAct_nonadapt(2) bestAct_nonadapt(3), bestAct_nonadapt(8) (bestAct_nonadapt(3)+[1:bestAct_nonadapt(4)]*bestAct_nonadapt(5)),...
            (bestAct_nonadapt(8)+[1:bestAct_nonadapt(9)]*bestAct_adapt(10))];
        
        % forward simulations of shortage and infrastructure costs
        totalCost_adapt = squeeze(totalCostTime_adapt(:,:,1:3))/1E6;% 1 is flex, 2 is static
        totalCost_nonadapt = squeeze(totalCostTime_nonadapt(:,:,1:3))/1E6;% 1 is flex, 2 is static
        damCost_adapt = squeeze(damCostTime_adapt(:,:,1:3))/1E6;% 1 is flex, 2 is static
        damCost_nonadapt = squeeze(damCostTime_nonadapt(:,:,1:3))/1E6;% 1 is flex, 2 is static
        
        for p = 1:length(P_ranges)
            % re-initialize mean cost arrays:
            totalCost_P = NaN(6,1);
            damCost_P = NaN(6,1);
            shortageCost_P = NaN(6,1);
            
            % calculate the average shortage, dam, and expansion costs for each
            % climate state:
            ind_P = ismember(P_state_adapt(:,5),P_ranges{p});
            
            % find climate subset from data for adaptive and non-adaptive
            % operations
            totalCost_P([3,1,5]) = squeeze(mean(sum(totalCost_nonadapt(ind_P,:,:),2)));
            damCost_P([3,1,5]) = squeeze(mean(damCost_nonadapt(ind_P,1,:)));
            shortageCost_P([2,1,3]) = squeeze(mean(sum(totalCost_nonadapt(ind_P,2:5,:)- damCost_nonadapt(ind_P,2:5,:),2)));
            
            totalCost_P([4,2,6]) = squeeze(mean(sum(totalCost_adapt(ind_P,:,:),2)));
            damCost_P([4,2,6]) = squeeze(mean(damCost_adapt(ind_P,1,:)));
            shortageCost_P([5,4,6]) = squeeze(mean(sum(totalCost_adapt(ind_P,2:5,:)- damCost_adapt(ind_P,2:5,:),2)));
            
            climate_cols = [[214, 84, 58];[110,110,110];[55,114,204]]/255; 
            
            hold on
            if cs == 1
            l(p,cs) = scatter([1,2,3,4.2,5.2,6.2],shortageCost_P,...
                100,climate_cols(p,:),'square','LineWidth',1.5);
            elseif cs == 2
                l(p,cs) = scatter([1,2,3,4.2,5.2,6.2],shortageCost_P,...
                100,climate_cols(p,:),'x','LineWidth',2);
            end
            
            hold on
            
            ax = gca;
            ax.LineWidth = 0.5;
            ax.YAxis.FontSize = 22;
        end
    end
    
end
    
    ylabel('Mean Shortage Cost (M$)\newline ','FontWeight','bold','FontSize',12)
    xlim([0.5,6.7])
    ylim([-8, 300])
    box on

    ax = gca;
    ax.XTick = [1,2,3,4.2,5.2,6.2];
    ax.XTickLabel = '';
    ax.LineWidth = 0.5;
    myLabels = { 'Static', 'Flexible', 'Flexible', 'Static', 'Flexible', 'Flexible'; 
        'Design', 'Design', 'Planning', 'Design', 'Design', 'Planning'};
    for i = 1:length(myLabels)
        text(ax.XTick(i), ax.YLim(1)-5, sprintf('%s\n%s', myLabels{:,i}), ...
            'horizontalalignment', 'center', 'verticalalignment', 'top');    
    end
    
    hold on
    mkrType(1) = scatter(-1,1,75,[0 0 0]/255,'square','LineWidth',1.5,'DisplayName',"c'_{alt}");
    hold on
    mkrType(2) = scatter(-1,1,75,[0 0 0]/255,'x','LineWidth',2,'DisplayName',"c'_{base}");
    
    leg1 = legend([mkrType(1) mkrType(2)],'Location','northeast','autoupdate','off','FontSize',10);
    title(leg1,"c' Scenario")
    set(leg1,'Box','off')
    
    hold on
    mkrColor(1) = bar(-1,1,'FaceColor',climate_cols(1,:),'DisplayName','Dry');
    hold on
    mkrColor(2) = bar(-1,1,'FaceColor',climate_cols(2,:),'DisplayName','Moderate');
    hold on
    mkrColor(3) = bar(-1,1,'FaceColor',climate_cols(3,:),'DisplayName','Wet');
    hold on
    ah1=axes('position',get(gca,'position'),'visible','off');
    
    leg2 = legend(ah1,[mkrColor(1) mkrColor(2) mkrColor(3)],'Location',...
        'northeast','autoupdate','off','FontSize',10,'LineWidth',1.5);
    title(leg2,'Final Climate')
    set(leg2,'Box','off','Position',leg2.Position + [0 -0.12 0 0])
    leg2.ItemTokenSize = [20,18];

    % set lower dual x-axis labels
    text(0.12, -0.11, "Static Operations",'FontSize',12,'FontWeight','bold')
    text(0.62, -0.11, "Flexible Operations",'FontSize',12,'FontWeight','bold')
    
    % panel label
    text(-0.17, 1.02, "(b)",'FontSize',12,'FontWeight','bold')

    % save figure as .pdf and .eps and update fontsize
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    cd(strcat(dir, 'Plots'))
    exportgraphics(gcf, 'Figures/mainbody/fig8_scatter.pdf', 'ContentType', 'vector'); 
    saveas(gcf, 'Figures/mainbody/fig8_scatter.eps')
    saveas(gcf, 'Figures/mainbody/fig8_scatter.jpg')




