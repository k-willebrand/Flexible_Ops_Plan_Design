%% RESERVOIR OPERATIONS SDP - SAVE SIMULATED MEAN 20-YEAR SHORTAGE COSTS

% DESCRIPTION:
% This is the main script utilized to optimize and simulate adaptive 
% (flexible) and non-adaptive (static) reservoir operations using SDP.

% This script is parallized. Multiple jobs can be submitted to high 
% performance computing clusters (e.g., across different reservoir design 
% capacities to effeciently optimize and simulate operating policies across
% different discrete climate states. Then, postprocessing is performed to 
% save results in a formate compatible with the infrastructure planning
% SDP code.

% SDP implemenation to design the optimal operating policy for different
% storage capacities for the Mwache Dam in Mombasa, Kenya.
% (adapted from M3O toolbox)

% STEPS:

% To run this script, update the reservoir capacity values you wish to 
% optimize for(can be an array), date, and filepath to save intermediate 
% results. Set the values for runParam, costParam, climParam, 
% and sys_param fields accordingly.

% After running this script, the post-processing script should be run prior
% to save shortage costs in a format compatible with the infrastructure
% planning SDP code.

% NOTE: For adaptive (flexible) operations, set runParam.adaptiveOps = true 
% else set runParam.adaptiveOps = false. Set runParam.adaptiveOps in
% outside and within the parallized loop.

%% Set file path

clear all;

if ~isempty(getenv('SLURM_JOB_ID'))
    projpath = '/home/users/keaniw/Fletcher_2019_Learning_Climate';
    jobid = getenv('SLURM_JOB_ID');
else
    % set local path
    projpath = 'C:/Users/keaniw/Documents/Research/Kenya Project/Project Code/Flexible_Ops_Plan_Design'; 
    jobid = 'na';
end

addpath(genpath(projpath))

% folder to save intermediate output files
date = '031122'; % set date for save name
fol = strcat('SDP_reservoir_ops/outputs'); 
mkdir(fol) % folder to save intermediate output files for future post-processing!

%% Set run parameters for shortage cost calculations

% Define reservoir capacities (can be an array of capacities)
% tip: submit multiple jobs with different paritions of entire capacity range 
storage_vals = [30:10:150]; % set reservoir capacities (MCM).

% descritized temperature and precipitation state space
s_T_abs = [26.25, 26.75, 27.25, 27.95, 28.8]; % deg. C
s_P_abs = 49:1:119; % [mm/month]; we later subset 66-97 mm/month

load('runoff_by_state_02Nov2021.mat'); % load T and P time series

% number of temperature and precipation states to calculate shortage costs
num_T_states = size(runoff,1); % temperature states
num_P_states = size(runoff,2); % precipitation states

if ~isempty(getenv('SLURM_JOB_ID'))
     poolobj = parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')));
     fprintf('Number of workers: %g\n', poolobj.NumWorkers)
end

% Define total reservoir capacity and dead storage (MCM)
for ss=1:length(storage_vals)
    
    storage = storage_vals(ss); % reservoir capacity (MCM)
    dead_storage = 20; % MCM
    
    for s_P = 1:num_P_states % set 18:49 to only consider 66-97 mm/mo
   
        P_state = s_P;
        
        %% calculate average shortage costs from the optimal policy for each climate state
        
        for s_T = 1:num_T_states
            
            T_state = s_T;
            disp(strcat('s: ', string(storage),' s_T:   ',string(s_T),' s_P:   ',string(s_P)));
            
            % Define adaptive or non-adaptive operations
            runParam = struct;
            sys_param = struct;
            grids = struct;
            policy = struct;
            climParam = struct;
            costParam = struct;
            
            % If true, then use climate adaptive inflow distributions in the SDP.
            % If false, use non-adaptive inflow distribution in the SDP.
            runParam.adaptiveOps = true;
            
            %% Set relevant sdp_climate.m main script parameters
            
            % Number of time periods
            runParam.N = 5; % Current SDP model requires N = 5
            
            % Number of years to generate in T, P, streamflow time series
            runParam.steplen = 200; 
            
            % Urban water demand scenarios (low = 150,000; high = 300,000)[m3/d](Fletcher 2019)
            runParam.domDemand = 186000; % 2020 design demand
            
            % If false, do not include deslination plant (default).
            runParam.desalOn = false;
            
            % If true, save results
            runParam.saveOn = true; % set to true when parellized
            
            % Number of T,P time series to generate using stochastic weather generator
            climParam.numSampTS = 21;
            
            % demand variables
            dmd_dom = cmpd2mcmpy(runParam.domDemand); % MCM/Y
            if runParam.desalOn
                dmd_dom = cmpd2mcmpy(300000); % high domestic demand scenario for desalination (MCM/Y)
            end
            dmd_ag = 12*[2.5 1.5 0.8 2.0 1.9 2.9 3.6 0.6 0.5 0.3 0.2 3.1]; % MCM/Y          
            demand = dmd_dom + dmd_ag; % MCM/Y
            
            % structure with system paramenters
            sys_param.simulation.delta = 1/12; % monthly time step (1/12 of a year)
            env_flow = 0; % MCM/Y
            sys_param.MEF = repmat(env_flow,12,1); % maximum environmental flow (MCM/Y)
            
            %% set SDP optimization
            
            % load synthetic inflow data
            if runParam.adaptiveOps % adaptive
                qq  = runoff{T_state,P_state,1}' ; % use current climate state data
            else % non-adaptive
                qq  = runoff{1,29,1}' ; % use initial climate state data
            end
            sys_param.simulation.q = qq ; % MCM/Y
            
            % discretization of variables' domain
            grids.discr_s = [0:1:(storage-dead_storage)]'; % MCM (effective storage)
            grids.discr_u = [[0:3:max(demand,[],'all')+env_flow+10],[125:100:3000],[3500:1000:9560]]' ; % MCM/Y
            grids.discr_q = [[0:3:max(prctile(runoff{1,71,1},80))],[max(prctile(runoff{1,71,1},80))+50:50:max(prctile(runoff{1,71,1},85))],[max(prctile(runoff{1,71,1},85))+50:1000:max(runoff{1,71,1},[],'all')+1000]]' ; % MCM/Y
            sys_param.algorithm = grids;
            sys_param.integration_substep = 300; % number of mass balance monthly substeps
            
            % set algorithm parameters
            sys_param.algorithm.name = 'sdp';
            sys_param.algorithm.Hend = 0 ; % penalty set to 0
            sys_param.algorithm.T = 12;    % the period is equal 1 year
            sys_param.algorithm.gamma = 1; % set future discount factor
            tol = -1;    % accuracy level for termination
            max_it = 10; % maximum iteration for termination
            
            % Estimate cyclostationary pdf (assuming log-normal distribution)
            T = sys_param.algorithm.T ;
            Ny = length(sys_param.simulation.q)/T*climParam.numSampTS;
            
            Q = reshape( sys_param.simulation.q, T, Ny ); % reshape inflow to monthly
            sys_param.algorithm.q_stat = nan(T,2) ;
            
            for i = 1:T
                qi = Q(i,:) ;
                sys_param.algorithm.q_stat(i,:) = lognfit(qi); % inflow distribution parameters for policy development
            end
            
            % create a table of expected monthly evaporation for each discrete storage state
            ee = NaN(length(grids.discr_s),T); % table containing expected evaporation rate by month and discrete storage state
            for e=1:length(grids.discr_s)
                if runParam.adaptiveOps % adaptive: use current climate state data
                    evap = evaporation_sdp(grids.discr_s(e)+dead_storage, T_ts{T_state,1},P_ts{P_state,1},climParam, runParam)';  % cyclostationary
                    ee(e,:) = mean(reshape( evap, T, Ny ),2)';
                else % non-adaptive: use initial climate state data
                    evap = evaporation_sdp(grids.discr_s(e)+dead_storage, T_ts{1,1},P_ts{29,1},climParam, runParam)';  % cyclostationary
                    ee(e,:) = mean(reshape( evap, T, Ny ),2)';
                end
            end
            
            sys_param.simulation.ev = ee ; % MCM/Y
            
            % compute minimum/maximum release for each climate state
            [vv, VV] = construct_rel_matrices(grids,sys_param); 
            sys_param.algorithm.min_rel = vv;
            sys_param.algorithm.max_rel = VV;
            
            %% run SDP optimization
            
            Hopt =  opt_sdp(tol, max_it, sys_param, runParam, costParam) ;
            policy.H = Hopt ;

            %% run simulation of SDP-policies
            
            % simulate policy using the actual climate state data
            qq  = runoff{T_state,P_state,1}' ; % 21 200-year simulations
            for e=1:length(grids.discr_s)
                evap = evaporation_sdp(grids.discr_s(e)+dead_storage, T_ts{T_state,1},P_ts{P_state,1},climParam, runParam)';  % cyclostationary
                ee(e,:) = mean(reshape( evap, T, Ny ),2)';
            end
            sys_param.simulation.q = qq ; % MCM/Y
            sys_param.simulation.ev = ee ; % MCM/Y
            
            % update steplen and numSampTS for simulation of policy
            runParam.steplen = 200; % years
            climParam.numSampTS = 21; % simulations
            
            q_sim = sys_param.simulation.q;
            s_init = 30; % set constant initial effective reservoir storage (MCM)
            J = NaN(runParam.steplen*T,climParam.numSampTS);
            s = NaN(runParam.steplen*T+1,climParam.numSampTS);
            r = NaN(runParam.steplen*T+1,climParam.numSampTS);
            ev = NaN(runParam.steplen*T+1,climParam.numSampTS); % evaporation
            dmd_unmet = NaN(runParam.steplen*T,climParam.numSampTS);
            dmd_unmet_ag = NaN(runParam.steplen*T,climParam.numSampTS);
            dmd_unmet_dom = NaN(runParam.steplen*T,climParam.numSampTS);
            
            disp('simulation running..')
            for i=1:climParam.numSampTS
                [J(:,i), s(:,i), r(:,i), ev(:,i), dmd_unmet(:,i), dmd_unmet_ag(:,i), dmd_unmet_dom(:,i)] = simulateGibe( q_sim(:,i), s_init, policy, sys_param, runParam, costParam);
            end
            
            objective_ts = J;
            storage_ts = s;
            release_ts = r;
            evaporation_ts = ev;
            unmet_ts = dmd_unmet;
            unmet_ag_ts = dmd_unmet_ag;
            unmet_dom_ts = dmd_unmet_dom;
            
            % calculate average 20-year steady state shortage cost (after first 20 years)
            SS_cumsum_obj = sum(objective_ts(T*20+1:end,:));
            num_SS_years = (length(objective_ts(:,1))-20*T)/(20*T); 
            shortageCost = mean(SS_cumsum_obj/num_SS_years);
            
            % save simulated result variables of interest
            if runParam.saveOn
                if runParam.adaptiveOps
                    filename = strcat(fol,'/adaptive_domagCost231_SSTest','_st',num2str(s_T),'_sp',num2str(s_P),'_s',string(storage),'_',date,'.mat');
                else
                    filename = strcat(fol,'/nonadaptive_domagCost231_SSTest','_st',num2str(s_T),'_sp',num2str(s_P),'_s',string(storage),'_',date,'.mat');
                end
                parsave(filename, shortageCost, objective_ts,  storage_ts, unmet_ag_ts, unmet_dom_ts);
            end
        end
    end
end