%% Calculatation of the climate transition matrix 

% Calculate the Bellman transition vector for the climate states using the
% Bayesian statistical model

%% Setup parameters

% Number of 20-year time periods
N = 5; 

% Set Emissions Scenario
emisScenario = {'RCP19' 'RCP26' 'RCP34' 'RCP45' 'RCP6' 'RCP7' 'RCP85'};
runParam.setPathway = emisScenario{7};

% Set up climate parameters
climParam = struct;

%  Number of simulations to use in order to estimate absolute T and P
%  values based on relative difference from one time period to the next
climParam.numSamp_delta2abs = 100000;

% Number of T,P time series to generate using stochastic weather generator
climParam.numSampTS = 21;

% Number of years to generate in T, P, streamflow time series
climParam.steplen = 100;

% If true, test number of simulated climate values are outside the range of
% the state space in order to ensure state space validity
climParam.checkBins = true;

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

%% Calculate and save the climate transition matrix 

load('BMA_results_RCP85_2020-11-14.mat')
[T_Precip, ~, ~] = bma2TransMat_P( NUP, s_P, N, climParam);
T_Temp = deterministicTempMatrix(N);  % temperature transitions modeled as deterministic
T_name = strcat('Data/T_Temp_Precip_', runParam.setPathway) % save a different transition matrix file for different emissions pathways

save(T_name, 'T_Temp', 'T_Precip')

