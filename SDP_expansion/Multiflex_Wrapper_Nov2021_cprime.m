
%% Wrapper Function to Find and Simulate Optimal Flexible Combinations

% This is a wrapper script for the infrastructure planning SDP.
% Refer to multiflex_sdp_climate_StaticFlex_DetT_Nov2021.m for parameter
% descriptions.

% For one set of transition matrices, do the following:
%   1. Find the optimal flexible design dam size
%   2. Find the optimal static dam size
%   4. Find the optimal flexible planned dam
%   3. Run the forward simulation with the optimal flex/static sizes


clear all;


if ~isempty(getenv('SLURM_JOB_ID'))
    projpath = '/home/users/keaniw/Fletcher_2019_Learning_Climate';
else
    projpath = '/Users/keaniw/Documents/Research/Kenya Project/Project Code/Flexible_Ops_Plan_Design';
end
addpath(genpath(projpath))

jobid = getenv('SLURM_JOB_ID');

% initialize cost model parameters (c' and discount rate can be an array)
cps = [1.25E-06]; 
discounts = [0];
x_6 = 5; %costParam.PercFlex: initial upfront capital cost increase flex design
x_7 = 1; %costParam.PercFlexExp: expansion cost of flexible dam  flex design (0)
x_11 = 0; %costParam.PercPlan: initial upfront capital cost increase (0);
x_12 = 50; %costParam.PercPlanExp: expansion cost of flexibly planned dam (0.5)

% preallocate SDP costs
val_static = zeros(11,2);
Vs_static = zeros(11,2);
Vd_static = zeros(11,2);
val_flex = zeros(10,2);
Vs_flex = zeros(10,2);
Vd_flex = zeros(10,2);
val_plan = zeros(10,2);
Vs_plan = zeros(10,2);
Vd_plan = zeros(10,2);

for d = 1:length(discounts)
    disc = discounts(d);
    
for c = 1:length(cps)
    c_prime = cps(c);

for z=1:2 % the adaptive and non-adaptive policies
    
    % 1. Find the optimal static dam size
    bestVal_static = inf;
    count = 1;
    
    for j=30:10:150
        
        stateMsg = strcat('z=',num2str(z),', j=', num2str(j));
        disp(stateMsg)
        
        x(1) = 2; %optParam.optFlex: (1)flexible design (2) static (3) flexible plan (0) policy;
        x(2) = j; % optParam.staticCap: static dam size [MCM]
        x(3) = 50; %optParam.smallFlexCap: unexpanded flexible design dam size [MCM]
        x(4) = 1;  %optParam.numFlex: number of possible expansion capacities [#]
        x(5) = 10; %optParam.flexIncr: increment of flexible expansion capacities [MCM]
        x(6) = x_6/100; %costParam.PercFlex: initial upfront capital cost increase flex design
        x(7) = x_7/100; %costParam.PercFlexExp: expansion cost of flexible dam  flex design (0)
        x(8) = 50; %optParam.smallPlanCap: unexpanded flexible plan dam size [MCM]
        x(9) = 1; %optParam.numPlan: number of possible expansion capacities [#]
        x(10) = 10; %optParam.planIncr:
        x(11) = x_11/100; %costParam.PercPlan: initial upfront capital cost increase (0);
        x(12) = x_12/100; %costParam.PercPlanExp: expansion cost of flexibly planned dam (0.5)
        x(13) = false; %runParam.forwardSim: (1) forward sim (0) no forward sim
        x(14) = z-1; %runParam.adaptiveOps: (1) adaptive operations (0) non-adaptive operations
        x(15) = disc/100; %0.03; %costParam.discountrate
        
        run('multiflex_sdp_climate_StaticFlex_DetT_Nov2021');
        val = V(1, 12, 1, 1)/1E6
        val_static(count,z) = V(1, 12, 1, 1)/1E6;
        Vs_static(count,z) = Vs(1, 12, 1, 1)/1E6;
        Vd_static(count,z) = Vd(1, 12, 1, 1)/1E6;
        count = count + 1;
        if val < bestVal_static
            bestVal_static = val;
            allV_static = V;
            bestAct_static = x;
            allX_static = X;
            allVs_static = Vs;
            allVd_static = Vd;
            
        end
    end
    
    bestAct_static(2)
    if z == 2
        file_name = strcat('Results/Results_SDP_expansion/BestStatic_adaptive_cp', regexprep(strrep(string(c_prime), '.', ''), {'-0'}, {''}),'_g7_percFlex',regexprep(strrep(string(x_6), '.', ''), {'-0'}, {''}),'_percExp',regexprep(strrep(string(x_7), '.', ''), {'-0'}, {''}),'_disc',string(disc),'_50PercExpCapOr150');
        %save(file_name, 'bestAct_static','bestVal_static','allV_static','allX_static','allVs_static','allVd_static');
    elseif z == 1
        file_name = strcat('Results/Results_SDP_expansion/BestStatic_nonadaptive_cp', regexprep(strrep(string(c_prime), '.', ''), {'-0'}, {''}),'_g7_percFlex',regexprep(strrep(string(x_6), '.', ''), {'-0'}, {''}),'_percExp',regexprep(strrep(string(x_7), '.', ''), {'-0'}, {''}),'_disc',string(disc),'_50PercExpCapOr150');
        %save(file_name, 'bestAct_static','bestVal_static','allV_static','allX_static','allVs_static','allVd_static');
    end
    
    
% 2. Find the optimal flexible design dam size
    bestVal_flex = inf;
    count = 1;
    
    for j=30:10:140
        for h=1 % flex expansion increments
            for g=1:7 % number of flex expansion(7-2*h)
                
                stateMsg = strcat('z=',num2str(z),', j=', num2str(j), ', g=', num2str(g), ', h=', ...
                    num2str(h));
                disp(stateMsg)
               
                if (j+h*10*g) > ceil(1.5*j/10)*10 % if false, then this combination is ok
                     disp('Flex design size not feasible');
                elseif (j+h*10*g) > 150  % place max on possible dam size
                    disp('Flex design size not feasible');
                elseif ((j+h*10*g) == ceil(1.5*j/10)*10) || ((j+h*10*g) == 150) % max initial dam 80 MCM
                    x(1) = 1; %optParam.optFlex: (1)flexible design (2) static (3) flexible plan (0) policy;
                    x(2) = 50; % optParam.staticCap: static dam size [MCM]
                    x(3) = j; %optParam.smallFlexCap: unexpanded flexible design dam size [MCM]
                    x(4) = g;  %optParam.numFlex: number of possible expansion capacities [#]
                    x(5) = h*10; %optParam.flexIncr: increment of flexible expansion capacities [MCM]
                    x(6) = x_6/100; %costParam.PercFlex: initial upfront capital cost increase flex design
                    x(7) = x_7/100; %costParam.PercFlexExp: expansion cost of flexible dam  flex design (0)
                    x(8) = 50; %optParam.smallPlanCap: unexpanded flexible plan dam size [MCM]
                    x(9) = 1; %optParam.numPlan: number of possible expansion capacities [#]
                    x(10) = 10; %optParam.planIncr:
                    x(11) = x_11/100; %costParam.PercPlan: initial upfront capital cost increase (0);
                    x(12) = x_12/100; %costParam.PercPlanExp: expansion cost of flexibly planned dam (0.5)
                    x(13) = false; %runParam.forwardSim: (1) forward sim (0) no forward sim
                    x(14) = z-1; %runParam.adaptiveOps: (1) adaptive operations (0) non-adaptive operations
                    x(15) = disc/100; %%costParam.discountrate
                    
                    run('multiflex_sdp_climate_StaticFlex_DetT_Nov2021');
                    val = V(1, 12, 1, 1)/1E6 % 1st temperature state, 77 mm/mo
                    val_flex(count,z) = V(1, 12, 1, 1)/1E6;
                    Vs_flex(count,z) = Vs(1, 12, 1, 1)/1E6;
                    Vd_flex(count,z) = Vd(1, 12, 1, 1)/1E6;
                    count = count + 1;
                    if val < bestVal_flex
                        bestVal_flex = val;
                        allV_flex = V;
                        bestAct_flex = x;
                        allX_flex = X;
                        allVs_flex = Vs;
                        allVd_flex = Vd;
                        
                    end
                else
                    disp('Flex design size not feasible');
                end
            end
        end
    end
    
    if z == 2
        file_name = strcat('Results/Results_SDP_expansion/BestFlex_adaptive_cp', regexprep(strrep(string(c_prime), '.', ''), {'-0'}, {''}),'_g7_percFlex',regexprep(strrep(string(x_6), '.', ''), {'-0'}, {''}),'_percExp',regexprep(strrep(string(x_7), '.', ''), {'-0'}, {''}),'_disc',string(disc),'_50PercExpCapOr150');
        %save(file_name, 'bestAct_flex','bestVal_flex','allV_flex','allX_flex','allVs_flex','allVd_flex');
    elseif z == 1
        file_name = strcat('Results/Results_SDP_expansion/BestFlex_nonadaptive_cp', regexprep(strrep(string(c_prime), '.', ''), {'-0'}, {''}),'_g7_percFlex',regexprep(strrep(string(x_6), '.', ''), {'-0'}, {''}),'_percExp',regexprep(strrep(string(x_7), '.', ''), {'-0'}, {''}),'_disc',string(disc),'_50PercExpCapOr150');
        %save(file_name, 'bestAct_flex','bestVal_flex','allV_flex','allX_flex','allVs_flex','allVd_flex');
    end
    
    
% 3. Find the optimal flexible planning dam size
    bestVal_plan = inf;
    count = 1;
    
    for j=30:10:140

        for h=1 % flex expansion increments
            for g=1:7 % number of flex expansion(7-2*h)

                stateMsg = strcat('z=',num2str(z),', j=', num2str(j), ', g=', num2str(g), ', h=', ...
                    num2str(h));
                disp(stateMsg)
                
                if (j+h*10*g) > ceil(1.5*j/10)*10 % if false, then this combination is ok
                     disp('Plan size not feasible');
                elseif (j+h*10*g) > 150 % place max on possible dam size
                    disp('Plan size not feasible');
                elseif ((j+h*10*g) == ceil(1.5*j/10)*10) || ((j+h*10*g) == 150) % if =50% initial capacity or 150 MCM
                    x(1) = 3; %optParam.optFlex: (1)flexible design (2) static (3) flexible plan (0) policy;
                    x(2) = 50; % optParam.staticCap: static dam size [MCM]
                    x(3) = 50; %optParam.smallFlexCap: unexpanded flexible design dam size [MCM]
                    x(4) = 1;  %optParam.numFlex: number of possible expansion capacities [#]
                    x(5) = 10; %optParam.flexIncr: increment of flexible expansion capacities [MCM]
                    x(6) = x_6/100; %costParam.PercFlex: initial upfront capital cost increase flex design
                    x(7) = x_7/100; %costParam.PercFlexExp: expansion cost of flexible dam  flex design (0)
                    x(8) = j; %optParam.smallPlanCap: unexpanded flexible plan dam size [MCM]
                    x(9) = g; %optParam.numPlan: number of possible expansion capacities [#]
                    x(10) = h*10; %optParam.planIncr:
                    x(11) = x_11/100; %costParam.PercPlan: initial upfront capital cost increase (0);
                    x(12) = x_12/100; %costParam.PercPlanExp: expansion cost of flexibly planned dam (0.5)
                    x(13) = false; %runParam.forwardSim: (1) forward sim (0) no forward sim
                    x(14) = z-1; %runParam.adaptiveOps: (1) adaptive operations (0) non-adaptive operations
                    x(15) = disc/100; %costParam.discountrate
                    
                    run('multiflex_sdp_climate_StaticFlex_DetT_Nov2021');
                    val = V(1, 12, 1, 1)/1E6
                    val_plan(count,z) = V(1, 12, 1, 1)/1E6;
                    Vs_plan(count,z) = Vs(1, 12, 1, 1)/1E6;
                    Vd_plan(count,z) = Vd(1, 12, 1, 1)/1E6;
                    count = count + 1;
                    if val < bestVal_plan
                        bestVal_plan = val;
                        allV_plan = V;
                        bestAct_plan = x;
                        allX_plan = X;
                        allVs_plan = Vs;
                        allVd_plan = Vd;
                        
                    end
                else
                    disp('Plan design size not feasible');
                end
            end
        end
    end

    if z == 2
        file_name = strcat('Results/Results_SDP_expansion/BestPlan_adaptive_cp', regexprep(strrep(string(c_prime), '.', ''), {'-0'}, {''}),'_g7_percFlex',regexprep(strrep(string(x_6), '.', ''), {'-0'}, {''}),'_percExp',regexprep(strrep(string(x_7), '.', ''), {'-0'}, {''}),'_disc',string(disc),'_50PercExpCapOr150');
        %save(file_name, 'bestAct_plan','bestVal_plan','allV_plan','allX_plan','allVs_plan','allVd_plan');
    elseif z == 1
        file_name = strcat('Results/Results_SDP_expansion/BestPlan_nonadaptive_cp', regexprep(strrep(string(c_prime), '.', ''), {'-0'}, {''}),'_g7_percFlex',regexprep(strrep(string(x_6), '.', ''), {'-0'}, {''}),'_percExp',regexprep(strrep(string(x_7), '.', ''), {'-0'}, {''}),'_disc',string(disc),'_50PercExpCapOr150');
        %save(file_name, 'bestAct_plan','bestVal_plan','allV_plan','allX_plan','allVs_plan','allVd_plan');
    end
        
% 4. Run the SDP and forward simulation with the optimal dam sizes
    
    % get optimal values
    x(1) = 0; %optParam.optFlex: (1)flexible design (2) static (3) flexible plan (0) policy;
    x(2) = bestAct_static(2); % optParam.staticCap: static dam size [MCM]
    x(3) = bestAct_flex(3); %optParam.smallFlexCap: unexpanded flexible design dam size [MCM]
    x(4) = bestAct_flex(4);  %optParam.numFlex: number of possible expansion capacities [#]
    x(5) = bestAct_flex(5); %optParam.flexIncr: increment of flexible expansion capacities [MCM]
    x(6) = x_6/100; %costParam.PercFlex: initial upfront capital cost increase flex design
    x(7) = x_7/100; %costParam.PercFlexExp: expansion cost of flexible dam  flex design (0)
    x(8) = bestAct_plan(8); %optParam.smallPlanCap: unexpanded flexible plan dam size [MCM]
    x(9) = bestAct_plan(9); %optParam.numPlan: number of possible expansion capacities [#]
    x(10) = bestAct_plan(10); %optParam.planIncr:
    x(11) = x_11/100; %costParam.PercPlan: initial upfront capital cost increase (0);
    x(12) = x_12/100; %costParam.PercPlanExp: expansion cost of flexibly planned dam (0.5)
    x(13) = true; %runParam.forwardSim: (1) forward sim (0) no forward sim
    x(14) = z-1; %runParam.adaptiveOps: (1) adaptive operations (0) non-adaptive operations
    x(15) = disc/100; %costParam.discountrate

    bestAct = x;
    
    run('multiflex_sdp_climate_StaticFlex_DetT_Nov2021');
    if z == 2
        file_name = strcat('Results/Results_SDP_expansion/BestFlexStaticPlan_adaptive_cp', regexprep(strrep(string(c_prime), '.', ''), {'-0'}, {''}),'_g7_percFlex',regexprep(strrep(string(x_6), '.', ''), {'-0'}, {''}),'_percExp',regexprep(strrep(string(x_7), '.', ''), {'-0'}, {''}),'_disc',string(disc),'_50PercExpCapOr150');
        %save(file_name, 'bestAct', 'V', 'X','action', 'totalCostTime','damCostTime','C_state','P_state','Vs','Vd')
    elseif z == 1
        file_name = strcat('Results/Results_SDP_expansion/BestFlexStaticPlan_nonadaptive_cp', regexprep(strrep(string(c_prime), '.', ''), {'-0'}, {''}),'_g7_percFlex',regexprep(strrep(string(x_6), '.', ''), {'-0'}, {''}),'_percExp',regexprep(strrep(string(x_7), '.', ''), {'-0'}, {''}),'_disc',string(disc),'_50PercExpCapOr150');
        %save(file_name, 'bestAct', 'V', 'X', 'action', 'totalCostTime','damCostTime','C_state','P_state','Vs','Vd')
    end
    
end
end
end