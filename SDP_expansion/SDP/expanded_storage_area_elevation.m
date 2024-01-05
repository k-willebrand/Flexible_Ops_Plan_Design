% DESCRIPTION: expanded area storage elevation by fitting a regression to the
% original model
% (compare results to the CES Feasibility study)

% load the storage-area-elevation data
load('area_storage.mat')
area = areaStorage(:,1);
storage = areaStorage(:,2);
elevation = [14, 15:5:90, 94, 95]';

% fit a line of best fit: elevation vs. capacity
fun = @(x, xdata) -x(1)*exp(-x(2)*xdata) + x(3);
%fun = @(x, xdata) x(1)*log(x(2)*xdata+x(3)) + x(4);
%p_sVe = polyfit(storage, elevation, 2);
x0 = [70, 0.01, 100];
x = lsqcurvefit(fun,x0,storage,elevation);

% Estimate elevation for larger dam sizes:
storage_new = [220:5:250]; % MCM
elev_new = fun(x, storage_new);

% confirm fit by plotting old data and new fitted points on plot
figure()
scatter(storage, elevation,'MarkerEdgeColor','b',...
    'DisplayName','Original Storage-Elevation Data')
hold on
plot([storage; storage_new'], fun(x, [storage; storage_new']),...
    'DisplayName',strcat("y= ",string(x(1)), " exp(", string(x(2)),"*x) + ", string(x(3))),'LineStyle','--','Color','r')
hold on
scatter(storage_new, elev_new,'Marker','*',...
    'DisplayName','Expanded Storage-Elevation Model','MarkerEdgeColor','r')
xlabel('Dam Storage (MCM)')
ylabel('Elevation (m)')
title('Dam Storage vs. Elevation for the Mwache Dam')
legend()
box on

%% Try fitting a cubic curve:
p = polyfit(storage, area, 3);

% Estimate elevation for larger dam sizes:
storage_new = [220:5:250]; % MCM
area_new = polyval(p, storage_new);

% confirm fit by plotting old data and new fitted points on plot
figure()
scatter(storage, area,'MarkerEdgeColor','b',...
    'DisplayName','Original Area-Storage Data')
hold on
plot([storage; storage_new'], polyval(p, [storage; storage_new']),...
    'DisplayName',strcat("y= ",string(p(1)), " x^3 + ", string(p(2)),"x^2 + ", string(p(3)), "x + ", string(p(4))),'LineStyle','--','Color','r')
hold on
scatter(storage_new, area_new,'Marker','*',...
    'DisplayName','Expanded Area-Storage Data','MarkerEdgeColor','r')
xlabel('Dam Storage (MCM)')
ylabel('Area (sq. km)')
title('Area vs. Dam Storage for the Mwache Dam')
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
costmodel = vertcat(costmodel, combine_new);
save('new_dam_cost_model.mat', 'costmodel')
