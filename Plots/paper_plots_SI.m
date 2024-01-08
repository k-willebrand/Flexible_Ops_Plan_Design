% Description: this script makes the figures used in the supplemental 
% information (SI) of the manuscript.

% Update to your local path
dir = '/Users/keaniw/Documents/Research/Kenya Project/Project Code/Flexible_Ops_Plan_Design/';
addpath(genpath(dir))

%% Figure S1: K-S Test for log-normal inflow assumption
% Description: Test that monthly inflow for each climate state follows a 
% log normal distribution that is statistically significant (K-S test)

% =========================================================================
% (1) LOAD THE DATA FILES TO PLOT

load('runoff_by_state_02Nov2021.mat');

% specificy climate states
s_T_abs = [26.25, 26.75, 27.25, 27.95, 28.8]; % deg. C
s_P_abs = 49:1:119; % mm/month ; 
s_P_ = 66:1:97;
ind_P = find(ismember(s_P_abs,s_P_) > 0); % we model 66 to 97 mm/month

% set algorithm parameters
T = 12;    % the period is equal 1 year or 12 months
numSampTS = 21; % number of simulations (21 GCMs)
steplen = 200; % number of years

% =========================================================================
% (2) CALCULATE P-VALUE FROM K-S TEST FOR EACH CLIMATE AND MONTH OF THE YEAR

h_table = cell(5, length(ind_P));
p_table = cell(5, length(ind_P));

h_table = cell(9, length(ind_P));
p_table = cell(9, length(ind_P));

% calculate K-S test statistics
for s_p=1:length(ind_P)
    
    P_state = ind_P(s_p);
    
    for s_t=1:5
        
        T_state = s_t;
        
        qq = runoff{T_state,P_state,1}';
        Ny = length(qq)/T*numSampTS;
   
        Q = reshape(qq, T, Ny);
        q_stat = nan(T,2);
        h = nan(T,1);
        p = nan(T,1);
        
        % perform K-S test for each month and climate state
        for i = 1:12 % for each month
            qi = Q(i,:);
            q_stat(i,:) = lognfit(qi);
            mu = q_stat(i,1);
            sigma = q_stat(i,2);
            r = lognrnd(mu,sigma);
            [h(i,1),p(i,1)] = kstest2(r,qi);
        end
        h_table{s_t, s_p} = h; % h parameter of ln inflow
        p_table{s_t,s_p} = p; % p-value of ln inflow
        
    end
end

% =========================================================================
% (3) BAR PLOT K-S TEST RESULTS (P-VALUE)(8 figures total)

s_ts = 1:5;

for i=4 % 8 figures
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    s_ps = [1+4*(i-1):4+4*(i-1)];
    s_ps = [1+8*(i-1):8+8*(i-1)]

    for p=1:length(s_ps)
        for t=1:length(s_ts)
            s_t = s_ts(t);
            s_p = s_ps(p);
            subplot(length(s_ps),length(s_ts),(p-1)*length(s_ts)+t);
            bar(1:12,p_table{s_t,s_p});
            hold on
            plot(0:13,ones(14)*0.05,'r', 'LineWidth',2);
            ylim([0,1]);
            if mod((p-1)*length(s_ts)+t,length(s_ts)) == 1
                label_p = ylabel({'P State: '; strcat(string(s_P_(s_p)),' mm/month')},'fontweight','bold','FontSize',12);
                label_p.Position(1) = -5; % change horizontal position of ylabel
                set(get(gca,'YLabel'),'Rotation',0)
            end
            if p == 1
                title(strcat('T State: ',string(s_T_abs(s_t)),' C'))
            end
            if p == length(s_ps) && i == 4
                xlabel('month of year', 'FontSize',12)
            end
            set(gca,'linewidth',2)
        end
    end
    %sgtitle({strcat('K-S Test p-values: Log Normal Inflow (', string(i), '/8)');'\alpha = 0.05'})

    % prepare figure to save
    origUnits = fig.Units;
    fig.Units = fig.PaperUnits; 
    fig.PaperSize = fig.Position(3:4);
    fig.Units = origUnits;

    % save figure as .pdf and .eps
    cd(strcat(dir, 'Plots'))
    exportgraphics(gcf, strcat('Figures/SI/figS1_',string(i),'_ksbar.pdf'))
    exportgraphics(gcf, strcat('Figures/SI/figS1_',string(i),'_ksbar.eps'))
    exportgraphics(gcf, strcat('Figures/SI/figS1_',string(i),'_ksbar.jpg'), 'Resolution', 800)
    exportgraphics(gcf, strcat('Figures/SI/figS1_',string(i),'_ksbar.tif'), 'Resolution', 800)

    % create single .pdf figure
    if i == 1
        exportgraphics(gcf, 'Figures/SI/figS1_ksbar.pdf')
    else
        exportgraphics(gcf,'Figures/SI/figS1_ksbar.pdf','Append',true)
    end

    close(fig)
end

% create single .pdf figure
fileNames = strcat('Figures/SI/figS1_',string(1:4),'_ksbar.tif');
f = figure;
montage(fileNames,"Size",[4 NaN], 'ThumbnailSize', [])
exportgraphics(gcf, 'Figures/SI/figS1_ksbar.pdf')

fileNames = strcat('figS1_',string(1:4),'_ksbar.tif');
out = imtile(strcat('figS1_',string(1:4),'_ksbar.tif'),'GridSize', [4, 1]);
f = figure;
imshow(out)
% prepare figure to save
origUnits = fig.Units;
fig.Units = fig.PaperUnits; 
fig.PaperSize = fig.Position(3:4);
fig.Units = origUnits;
exportgraphics(gcf, 'Figures/SI/figS1_ksbar.pdf', 'Resolution', 1200)

%% Figure S2: Box plot of cost savings for flex design and planning by climate

% Description: Comparison of cost savings for flexible design vs. flexible 
% planning. The first column shows the frequency that flexible design 
% performs better than flexible planning in minimizing total lifetime costs
% and vice versa across all 10,000 Monte Carlo simulations, and subsets of 
% the 10,000 Monte Carlo simulations for representative dry, moderate, and 
% wet final climates assuming flexible operations. The second column shows 
% a box plot of simulated total cost savings between flexible design and 
% flexible planning in minimizing total cost and vice versa under flexible 
% operations.

% =========================================================================

% initialize the final figure figure
f=figure('Position', [476,329,560,451]);
t = tiledlayout(4,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

facecolors = [[153,204,204]; [187,221,221]; [153, 153, 204]; [204, 204, 255];...
    [255, 102, 102]; [255, 153, 153]]/255;

% =========================================================================
% (1) LOAD THE DATA FILES TO PLOT

cd('/Users/keaniw/Documents/Research/Kenya Project/Project Code/Flexible_Ops_Plan_Design/')

% specify file path for data
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

% make a generalized vector of 10,000 indices
allIdx = 1:length(totalCostTime_adapt);

% consider only flexible operations indices for flex planning and design
TotalCostDiff_adapt = (sum(totalCostTime_adapt(:,:,1),2) ...
    - sum(totalCostTime_adapt(:,:,3),2))/1E6; % cost flex design - flex plan
Idx_bestPlan_adapt = allIdx(TotalCostDiff_adapt > 0); % flex plan better
Idx_bestFlex_adapt = allIdx(TotalCostDiff_adapt < 0); % flex design better

% initialize the final figure figure
t = tiledlayout(4,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

facecolors = [[153,204,204]; [187,221,221]; [153, 153, 204]; [204, 204, 255];...
    [255, 102, 102]; [255, 153, 153]]/255;

% =========================================================================
% (2) PLOT COST SAVINGS FOR DIFFERENT CLIMATE SUBSETS OF THE SIMULATIONS

P_ranges = {[66:97], [66:76],[77:86],[87:97]}; % all, dry, moderate, and wet climates

for p = 1:length(P_ranges)
    
    % index climate scenario
    ind_P = ismember(P_state_adapt(:,5),P_ranges{p});

    % consider only flexible operations indices for flex planning and design
    TotalCostDiff_adapt = (sum(totalCostTime_adapt(ind_P,:,1),2) ...
        - sum(totalCostTime_adapt(ind_P,:,3),2))/1E6; % cost flex design - flex plan
    Idx_bestPlan_adapt = allIdx(TotalCostDiff_adapt > 0); % flex plan better
    Idx_bestFlex_adapt = allIdx(TotalCostDiff_adapt < 0); % flex design better

    
    % (a) BAR GRAPH: NUMBER SIMS FLEX DESIGN VS. FLEX PLANNING PERFORM BETTER

    % flexible ops
    nexttile([1,1])
    b = bar([length(Idx_bestFlex_adapt); length(Idx_bestPlan_adapt)],'FaceColor','flat');
    for i=1:length(b.XData)
        b.CData(i,:) = facecolors(2+2*i,:);
    end
    xticklabels({'Flexible Design', 'Flexible Planning'})
    set(gca, 'YGrid', 'on', 'XGrid', 'off');
    xlim([0.5,2.5])
    ylim([0,9000])
    if p == 4
    else
        set(gca,'XTick',[])
        if p == 1
            title('Frequency of Best Performance\newline ','FontWeight','bold', 'FontSize', 13)
        end
    end

    % (b) BOXPLOT: TOTAL COST SAVINGS FROM SIMULATIONS
    
    % specify file naming convensions
    disc = string(0);
    cp = 1.25e-6;
    c_prime = regexprep(strrep(string(cp), '.', ''), {'-0'}, {''});
             
    % load flexible operations files
    load(strcat(folder,'/BestFlexStaticPlan_adaptive_cp',c_prime,'_g7_percFlex5_percExp1_disc',disc,'_50PercExpCapOr150.mat'))
    
    % specify variable names
    totalCostTime_adapt = totalCostTime;
    P_state_adapt = P_state;
    
    %P_ranges = {[66:76],[77:86],[87:97]}; % reference dry, moderate, and wet climates
    
    % make a generalized vector of 10,000 indices
    allIdx = 1:length(totalCostTime_adapt);
    
    % flexible operations
    nexttile([1,1])

    if isempty(TotalCostDiff_adapt(Idx_bestFlex_adapt)) == 0 && isempty(TotalCostDiff_adapt(Idx_bestPlan_adapt)) == 0
        groups = [repmat({'Flexible Design'},length(Idx_bestFlex_adapt),1); ...
            repmat({'Flexible Planning'},length(Idx_bestPlan_adapt),1)];
        boxplot([-TotalCostDiff_adapt(Idx_bestFlex_adapt);TotalCostDiff_adapt(Idx_bestPlan_adapt)],...
            groups,'BoxStyle','filled','Widths',0.3,'OutlierSize',5,'Symbol','.','BoxStyle','filled')
    elseif isempty(TotalCostDiff_adapt(Idx_bestFlex_adapt))
        groups = [repmat({'Flexible Design'},1,1); ...
            repmat({'Flexible Planning'},length(Idx_bestPlan_adapt),1)];
        boxplot([0;TotalCostDiff_adapt(Idx_bestPlan_adapt)],...
            groups,'BoxStyle','filled','Widths',0.3,'OutlierSize',5,'Symbol','.','BoxStyle','filled')
    elseif isempty(TotalCostDiff_adapt(Idx_bestPlan_adapt))
        groups = [repmat({'Flexible Design'},length(Idx_bestFlex_adapt),1); ...
            repmat({'Flexible Planning'},1,1)];
        boxplot([-TotalCostDiff_adapt(Idx_bestFlex_adapt);0],...
            groups,'BoxStyle','filled','Widths',0.3,'OutlierSize',5,'Symbol','.','BoxStyle','filled')
    end
    hold on

    plot(1,mean(-TotalCostDiff_adapt(Idx_bestFlex_adapt)), 'k.', 'MarkerSize',15)
    plot(2,mean(-TotalCostDiff_adapt(Idx_bestPlan_adapt)), 'k.', 'MarkerSize',15)
    whisks = findobj(gca,'Tag','Whisker');
    outs = findobj(gca, 'type', 'line','Tag', 'Outliers');
    meds = findobj(gca, 'type', 'line', 'Tag', 'Median');
    set(meds(1),'Color','k', 'LineWidth', 1.25);
    set(meds(2),'Color','k', 'LineWidth', 1.25);
    set(whisks(1),'Color',facecolors(6,:), 'LineWidth', 2);
    set(whisks(2),'Color',facecolors(4,:), 'LineWidth', 2);
    set(outs(1),'MarkerEdgeColor',facecolors(6,:), 'MarkerSize', 10);
    set(outs(2),'MarkerEdgeColor',facecolors(4,:),'MarkerSize',10);
    a = findobj(gca,'Tag','Box');
    set(a(1),'Color',facecolors(6,:),'Linewidth',25);
    set(a(2),'Color',facecolors(4,:),'Linewidth',25);
    ylim([0,17])
    set(gca, 'YGrid', 'on', 'XGrid', 'off', 'YAxisLocation','right');
    title(' \newline ','FontWeight','bold')

    if p == 4
    else
        set(gca,'XTick',[])
        if p == 1
            title('Simulated Total Cost Savings\newline ','FontWeight','bold', 'FontSize', 13)
        end
    end
end

% add climate titles
text(-0.25, 104, 'All Final Climates', 'FontSize',13, 'FontWeight', 'bold');
text(-0.25, 76, 'Dry Final Climates', 'FontSize',13, 'FontWeight', 'bold');
text(-0.44, 48, 'Moderate Final Climates', 'FontSize',13, 'FontWeight', 'bold');
text(-0.25, 20, 'Wet Final Climates', 'FontSize',13, 'FontWeight', 'bold');

% add axis titles
ylabel(t, 'Number of Simulations','FontWeight','bold');
z = text(2.7, 36, 'Cost Savings ($M)', 'FontWeight','bold', 'FontSize', 12);
z.Rotation = 90;

% save figure as .pdf and .eps
cd(strcat(dir, 'Plots'))
exportgraphics(gcf, 'Figures/SI/figS2_barbox.pdf', 'ContentType', 'vector'); 
saveas(gcf, 'Figures/SI/figS2_barbox.eps')
saveas(gcf, 'Figures/SI/figS2_barbox.jpg')


