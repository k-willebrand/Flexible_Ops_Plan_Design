function [ H , idx_u ] = Bellman_sdp( H_ , s_curr, moy,  sys_param, runParam, costParam)

%global sys_param runParam costParam

%-- Initialization --
discr_s = sys_param.algorithm.discr_s;
discr_u = sys_param.algorithm.discr_u;
discr_q = sys_param.algorithm.discr_q;
n_u = length(sys_param.algorithm.discr_u);
n_q = length(sys_param.algorithm.discr_q);
gamma   = sys_param.algorithm.gamma;
delta = sys_param.simulation.delta;

VV = repmat(sys_param.simulation.VV', 1 , n_u);
vv = repmat(sys_param.simulation.vv', 1 , n_u);

mi_q    = sys_param.algorithm.q_stat(moy,1); % mean of log(disturbance)
sigma_q = sys_param.algorithm.q_stat(moy,2); % std of log(disturbance)

%-- Calculate actual release contrained by min/max release rate --
uu = repmat(discr_u', n_q, 1);
R = min( VV , max( vv , uu ) ); % possible releases constrained by min and max releases

%==========================================================================
% Calculate the state transition; TO BE ADAPTAED ACCORDING TO
% YOUR OWN CASE STUDY
evap = interp1(discr_s, sys_param.simulation.ev(:,moy),s_curr,'pchip',max(sys_param.simulation.ev(:,moy))); % evap in MCM/Y: lookup expected evap based on incremental storage
%evap = interp1(discr_s, sys_param.simulation.ev(:,moy),s_curr,'linear','extrap');

qq = repmat( discr_q, 1, n_u ); % MCM/Y
s_eff = max(s_curr + delta*qq - delta*evap,0); % effective storage
%s_eff = max(s_curr + delta*qq - delta*evap(moy),0); % effective storage
s_next = s_eff - delta * R; % MCM


%==========================================================================
% Compute immediate costs and aggregated one; TO BE ADAPTAED ACCORDING TO
% YOUR OWN CASE STUDY

% Define target demands
dmd_dom = cmpd2mcmpy(runParam.domDemand); % MCM/Y
if runParam.desalOn
    dmd_dom = cmpd2mcmpy(300000); % high domestic demand scenario for desal
end
dmd_ag = 12*[2.5 1.5 0.8 2.0 1.9 2.9 3.6 0.6 0.5 0.3 0.2 3.1]; % MCM/Y (Fletcher)
%dmd_ag = [32.65 35.40 30.31 10.47 3.44 9.09 13.78 19.70 23.42 21.22 21.22 24.11]; % MCM/Y (MWI Average rescaled)
demand = dmd_dom + dmd_ag; % MCM/Y

% Calclate unmet demands from releases and this period's shortage cost (G)
% ag demands are unmet before dom demands
unmet = max((demand(moy) - R)*delta*1E6, 0); % CM
unmet_ag = min(unmet, dmd_ag(moy)*delta*1E6); % CM
unmet_dom = unmet - unmet_ag; % CM
%G =  (unmet_ag.^2 * costParam.agShortage + unmet_dom.^2 * costParam.domShortage); % CUBIC METERS ^ 2 forumulation
G =  (unmet_ag).^2 + 2.3.*(unmet_dom).^2; % CUBIC METERS ^ 2 forumulation

%-- Compute cost-to-go given by Bellman function --
H_ = interp1( discr_s , H_ , s_next(:), 'linear', 'extrap' ) ;
H_ = reshape( H_, n_q, n_u );

%-- Compute resolution of Bellman value function --
% compute the probability of occurence of inflow that falls within the
% each bin of descritized inflow level
cdf_q      = logncdf( discr_q , mi_q , sigma_q );  
p_q        = diff(cdf_q);                          
p_diff_ini = 1-sum(p_q);                           
p_diff     = [ p_diff_ini ; p_q];
cdf_diff = cumsum(p_diff); % cdf for discretized inflow

% for expected value:
Q     = (G + gamma.*H_)'*p_diff ; % expected value
H     = min(Q)                  ;
sens  = eps                     ;
idx_u = find( Q <= H + sens )   ;

% for percentile optimization:
prctile_opt = 0.1; % percentile to optimize for
Q_ini = (G + gamma.*H_) ; 
% value at desired inflow percentile
Q = Q_ini(abs(cdf_diff-prctile_opt) == min(abs(cdf_diff-prctile_opt)),:);
H = min(Q);
sens  = eps;
idx_u = find(Q <= H + sens);

% Q_max = max(G + gamma.*H_); % robust, worst-case value
% Q_min = min(G + gamma.*H_); % optimistic, best-case value
% Q_prct = prctile(G + gamma.*H_,90)'; % 90th percentile
% H     = min(Q_prct)             ;
% diffValues = Q_prct - H;
% diffValues(diffValues < 0) = inf;
% idx_u = find(Q_prct == min(diffValues));

end
