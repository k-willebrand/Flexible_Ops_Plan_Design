% DESCRIPTION: expanded dam cost model by fitting a regression to the
% original dam cost model between storage and total dam cost

% load the data
load('dam_cost_model.mat')

% fit a line of best fit
p = polyfit(costmodel.storage, costmodel.dam_cost, 2);

% apply the line of best fit for larger dam sizes:
storage_new = [20:5:45]; % MCM
dcost_new = polyval(p, storage_new);

% confirm fit by plotting old data and new fitted points on plot
figure()
scatter(costmodel.storage, costmodel.dam_cost,'MarkerEdgeColor','b',...
    'DisplayName','Original Dam Cost Model Data')
hold on
scatter(storage_new, dcost_new,'Marker','*',...
    'DisplayName','Expanded Dam Cost Model Data','MarkerEdgeColor','r')
hold on
plot([storage_new'; costmodel.storage], [dcost_new'; costmodel.dam_cost],...
    'DisplayName',strcat("y=",string(p(1)),"x^2 + ", string(p(2)),"x + ",...
    string(p(3))),'LineStyle','--','Color','r')
xlabel('Dam Storage (MCM)')
ylabel('Dam Cost ($)')
title('Dam Cost vs. Storage for the Mwache Dam')
legend()
box on

%% Create a new dam cost model file for the expanded dam sizes

% Calculate the additional dam cost model parameters:
unit_cost_new = dcost_new./storage_new;
exp_unit_cost_new = 1.7*unit_cost_new; % based on the original file's convention
height_new = nan(1, length(unit_cost_new));

% Create a table of the new cost model paramters:
combine_new = array2table([storage_new', height_new', dcost_new', unit_cost_new',...
    exp_unit_cost_new'],'VariableNames',costmodel.Properties.VariableNames);

% combine new dam costs with old dam costs and save new file
costmodel = vertcat(combine_new, costmodel);
save('new_dam_cost_model.mat', 'costmodel')
