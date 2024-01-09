%% Calculatation of the temperature, precipitation, and runoff monthly time series

% Use k-nn stochastic weather generator (Rajagopalan et al. 1999) to
% generate time series of monthly T and P based on 20-year means from state
% space.  Then, use CLIRUN hydrological model to simulate runoff monthly time 
% series for each T,P time series


%% Setup 

% Set Project root folder and Add subfolders to path; runs either on desktop 
% or on a cluster using SLURM queueing system 
if ~isempty(getenv('SLURM_JOB_ID'))
    projpath = '/net/fs02/d2/sfletch/Mombasa_climate';
else
    projpath = '/Users/sarahfletcher/Dropbox (MIT)/Fletcher_2019_Learning_Climate';
end
addpath(genpath(projpath))

jobid = getenv('SLURM_JOB_ID');

% Get date for file name when saving results 
datetime=datestr(now);
datetime=strrep(datetime,':','_'); %Replace colon with underscore
datetime=strrep(datetime,'-','_');%Replace minus sign with underscore
datetime=strrep(datetime,' ','_');%Replace space with underscore

%% Setup parameters

% Number of 20-year time periods
N = 5; 

runParam = struct;

% Set Emissions Scenario
emisScenario = {'RCP19' 'RCP26' 'RCP34' 'RCP45' 'RCP6' 'RCP7' 'RCP85'};
runParam.setPathway = emisScenario{7};

% If true, simulate runoff time series from T, P time series using CLIRUN. If false, load saved.
runParam.runRunoff = true; 

% If true, simulate T, P time series from mean T, P states using stochastic weather gen. If false, load saved.
runParam.runTPts = true; 

% If true, change indices of saved runoff time series to correspond to T, P states (needed for parfor implementation)
runParam.runoffPostProcess = true; 

% Set up climate parameters
climParam = struct;

%  Number of simulations to use in order to estimate absolute T and P
%  values based on relative difference from one time period to the next
climParam.numSamp_delta2abs = 100000;

% Number of T,P time series to generate using stochastic weather generator
climParam.numSampTS = 21;

% Number of years to generate in T, P, streamflow time series
climParam.steplen = 100;

% Define state space for mean 20-year precipitation and temperature
% output: 66 mm/month to 97 mm/month

% Percent change in precip from one time period to next
climParam.P_min = -.3;
climParam.P_max = .3;
climParam.P_delta = .02; 
s_P = climParam.P_min : climParam.P_delta : climParam.P_max;
climParam.P0 = s_P(15);
climParam.P0_abs = 77; %mm/month
M_P = length(s_P);

% Change in temperature from one time period to next
% output: 26.25, 26.75, 27.25, 27.95, 28.8 deg. C

climParam.T_delta = [0.5 0.5 0.7 0.85];
climParam.T0_abs = 26.25;
M_T = N;

% Absolute temperature values
T_abs_max = sum(climParam.T_delta);
s_T_abs(1) = climParam.T0_abs;
for i=1:N-1
    s_T_abs(i+1) = s_T_abs(i) + climParam.T_delta(i);
end
M_T_abs = length(s_T_abs);
T_Temp_abs = zeros(M_T_abs,M_T_abs,N);

% Absolute precip values
P_abs_max = max(s_P) * N;
s_P_abs = 66:1:97;
M_P_abs = length(s_P_abs);
P_bins = [s_P_abs-climParam.P_delta/2 s_P_abs(end)+climParam.P_delta/2];
T_Precip_abs = zeros(M_P_abs,M_P_abs,N);



%% Load transition probability matrix

% load pre-saved transition probability matrix
load('T_Temp_Precip_RCP85B', 'T_Precip')
T_Temp = deterministicTempMatrix(N);

% Prune state space -- no need to calculate policies for T and P states
% that are never reached when simulating future climates based on Bayesian
% model
for t = 1:N
    index_s_p_time{t} = find(~isnan(T_Precip(1,:,t)));
    index_s_t_time{t} = find(~isnan(T_Temp(1,:,t)));
end

%% Calculate and save T, P monthly time series for each long-term T, P state

% Use k-nn stochastic weather generator (Rajagopalan et al. 1999) to
% generate time series of monthly T and P based on 20-year means from state
% space

T_ts = cell(M_T_abs,N);
P_ts = cell(M_P_abs,N);

[Tanom, Panom] = mean2TPtimeseriesMJL_2(1, climParam.steplen, climParam.numSampTS); 
for t = 1:N

    for i = 1:M_T_abs  
        T_ts{i,t} = Tanom + s_T_abs(i)*ones(size(Tanom));
    end

    for i = 1:M_P_abs  
        Ptmp = Panom + s_P_abs(i)*ones(size(Tanom));
        Ptmp(Ptmp<0) = 0;
        P_ts{i,t} = Ptmp;
    end

end

savename_runoff = strcat('Data/runoff_by_state_', jobid,'_', datetime);
save(savename_runoff, 'T_ts', 'P_ts')

%% Calculate and save runoff monthly time series for each long-term T, P state

% Use CLIRUN hydrological model to simulate runoff monthly time series for
% each T,P time series

if runParam.runRunoff 

    % Generate runoff timeseries - different set for each T,P combination
    runoff = cell(M_T_abs, M_P_abs, N);


    % Set up parallel for running on cluster with SLURM queueing system
    if ~isempty(getenv('SLURM_JOB_ID'))
        poolobj = parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')));
        fprintf('Number of workers: %g\n', poolobj.NumWorkers)
    end

    for t = 1

        % loop over available temp states
        index_s_t_thisPeriod = index_s_t_time{t}; 
        parfor i = 1:length(index_s_t_thisPeriod)
            index_s_t = index_s_t_thisPeriod(i);

            runoff_temp = cell(M_P_abs,1);

            % loop over available precip states
            index_s_p_thisPeriod = index_s_p_time{t}; 
            for index_s_p = index_s_p_thisPeriod

                disp(strcat('i: ', string(i), '   index_s_p: ', string(index_s_p)))
                % Call CLIRUN streamflow simulator
                runoff_temp{index_s_p} = ...
                    TP2runoff(T_ts{index_s_t,t}, P_ts{index_s_p,t}, climParam.steplen);

            end

            runoff(i, :, t) = runoff_temp;

        end
    end


    savename_runoff = strcat('Data/runoff_by_state_', jobid,'_', datetime);
    save(savename_runoff, 'runoff', 'T_ts', 'P_ts')


    if runParam.runoffPostProcess
        % The nature of the parfor loop above saves the runoff timeseries in
        % first available index; this section moves to correct cell
        % corresponsing to P, T state space
        runoff_post = cell(M_T_abs, M_P_abs, N);
        for t = 1:N

            index_s_p_thisPeriod = index_s_p_time{t}; 
            for index_s_p = index_s_p_thisPeriod

                index_s_t_thisPeriod = index_s_t_time{t}; 
                for i= 1:length(index_s_t_thisPeriod)

                    runoff_post{index_s_t_thisPeriod(i),index_s_p,t} = runoff{i,index_s_p,t};

                end
            end
        end

        runoff = runoff_post;

        savename_runoff = strcat('Data/runoff_by_state_', jobid,'_', datetime);
        save(savename_runoff, 'runoff', 'T_ts', 'P_ts')

    end
end