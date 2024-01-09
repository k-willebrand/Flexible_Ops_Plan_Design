
function [T_Precip_abs, T_Precip_delta, P_over_time] ...
    = bma2TransMat_P(NUP, s_P, N, climParam)

% Inputs BMA samples and returns Bellman transition matrices

P_abs_max = max(s_P) * N;
s_P_abs = 66:1:97;
M_P_abs = length(s_P_abs);
P_bins_abs = [s_P_abs-climParam.P_delta/2  s_P_abs(end)+climParam.P_delta/2];
T_Precip_abs = zeros(M_P_abs,M_P_abs,N);


%% Deltas

P_step = s_P(2) - s_P(1);
P_bins = [s_P-P_step/2 s_P(end)+P_step/2];
M_P = length(s_P);

% Check how many samples outside bins
if climParam.checkBins

    numPSamp = numel(NUP);
    indP = find(NUP < P_bins(1) | NUP > P_bins(end));
    percPout = length(indP)/numPSamp
    
end
%% If there are values outside of the range, round them

% Round NUT values outside of the range into the range

% Round NUP values outside of the range into the range
if percPout > 0
    for i=1:1000
        for t = 1:N
            for index_s_P = 1:M_P
                if NUP(i,t,index_s_P) < P_bins(1);
                    NUP(:,t,index_s_P) = P_bins(1);
                elseif NUP(i,t,index_s_P) > P_bins(end)
                    NUP(i,t,index_s_P) = P_bins(end);
                end 
            end
        end
    end
end

% Check that issue is fixed

numPSamp = numel(NUP);
indP = find(NUP < P_bins(1) | NUP > P_bins(end));
percPout2 = length(indP)/numPSamp

%% Calculate transition matrices for deltas

T_Precip_delta = zeros(M_P, M_P, N);
for t = 1:N
    for index_s_P = 1:M_P
        T_Precip_delta(:,index_s_P,t) =  histcounts(NUP(:,t,index_s_P), P_bins, 'Normalization', 'Probability');
    end
end

%% Time series from deltas

% Simulate delta time seriess
p = rand(climParam.numSamp_delta2abs,N+1);
state_ind_P = zeros(climParam.numSamp_delta2abs,N);
P0_ind = randi(M_P,climParam.numSamp_delta2abs,1);

state_ind_P(:,1) = P0_ind;
for i = 1:climParam.numSamp_delta2abs
    for t = 1:N
        state_ind_P(i,t+1) = find(p(i,t) < cumsum(T_Precip_delta(:,state_ind_P(i,t),t)),1);
    end
end
P_delta_over_time = s_P(state_ind_P);


% Select starting point 
P0_abs_ind = randi(M_P_abs,climParam.numSamp_delta2abs,1); %I think we'd still want rounding for precip
P0_abs = s_P_abs(P0_abs_ind)';

% Precip is percent change
P_over_time = zeros(climParam.numSamp_delta2abs, N+1);

for i = 1:climParam.numSamp_delta2abs
    for t = 1:6
        if t == 1 
            P_over_time(i,t) = P0_abs(i);
        else
            P_over_time(i,t) = P_over_time(i,t-1) * (1+P_delta_over_time(i,t));
        end
    end
end


%% Round values so that temp and precip values do not extend past state space

% Check number of values initially outside
numPSamp = numel(P_over_time);
indP = find(P_over_time < s_P_abs(1) | P_over_time > s_P_abs(end));
percPout = length(indP)/numPSamp

%Round P_over_time values outside of range
if percPout > 0
    for i=1:100000
        for t = 1:N+1
            if P_over_time(i,t) < s_P_abs(1);
                P_over_time(i,t) = s_P_abs(1);
            elseif P_over_time(i,t) > s_P_abs(end)
                P_over_time(i,t) = s_P_abs(end);
            end
        end
    end
end

% Check that issue is fixed
numPSamp = numel(P_over_time);
indP = find(P_over_time < s_P_abs(1) | P_over_time > s_P_abs(end));
percPout2 = length(indP)/numPSamp

%% Absolutes from time series

% Calculate conditional prob of going from state X in time t-1 to state Y
% in time t

P_over_time_rounded = round(P_over_time);
for i = 1:length(s_P_abs)
    for t = 1:N
        P_current = s_P_abs(i);
        indexNow = find(P_over_time_rounded(:,t) == P_current);
        relevant_deltas = P_over_time_rounded(indexNow,t+1) - P_over_time_rounded(indexNow,t);
        P_next = P_current + relevant_deltas;
        T_Precip_abs(:,i,t) = histcounts(P_next, P_bins_abs, 'Normalization', 'Probability');     
    end
end

%% Get rid of NaN

temp_T_Precip = T_Precip_abs;
ind = find(isnan(T_Precip_abs));
temp_T_Precip(ind) = 0;
[ind1, ind2, ind3] = ind2sub(size(sum(temp_T_Precip)), find(sum(temp_T_Precip) == 0));
index = sub2ind(size(temp_T_Precip), ind2, ind2, ind3);
temp_T_Precip(index) = 1; %stay in this state
if sum(find( abs(sum(temp_T_Precip) -1) > .001)) > 0
    error('invalid T Precip')
end
T_Precip_abs = temp_T_Precip;


