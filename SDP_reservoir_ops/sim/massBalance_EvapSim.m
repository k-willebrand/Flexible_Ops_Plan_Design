function [s1,r1] = massBalance_EvapSim( s, u, q, mo ,discr_s, sys_param, T, P)
% Output:
%      s1 - final storage. 
%      r1 - release over the 24 hours. 
%
% Input:
%       s - initial storage. 
%       u - release decision 
%       q - inflow 
%       mo - cumulative month of time series
%
% See also MIN_RELEASE, MAX_RELEASE

%global sys_param

HH = sys_param.integration_substep;
delta = sys_param.simulation.delta/HH; % Y
s_ = nan(HH+1,1);
r_ = nan(HH+1,1);
evap = interp1(discr_s, sys_param.simulation.ev(:,moy),s,'pchip',max(sys_param.simulation.ev(:,moy))); % evap in MCM/Y: lookup expected evap based on initial storage
%evap = sys_param.simulation.ev  ; % evap in MCM/Y

s_(1) = s;
s_max = discr_s(end); % maximum effective reservoir storage
ee = evaporation_sdp(s,T,P,climParam,runParam);
for i=1:HH
    %evap = interp1(discr_s, sys_param.simulation.ev(:,moy),s_(i),'pchip',max(sys_param.simulation.ev(:,moy))); % evap in MCM/Y: lookup expected evap based on incremental storage
    %s_eff = max(s_(i) + delta*q - delta*evap(moy),0); % effective storage
    s_eff = max(s_(i) + delta*q - delta*evap,0); % effective storage
    [r_min, r_max] = min_max_relG3(s_eff, moy, s_max, sys_param) ;
    r_(i+1) = min( r_max , max( r_min , u ) ); % MCM/Y
    s_(i+1) = s_eff - delta*r_(i+1); % MCM
end

s1 = s_(HH);
r1 = mean(r_(2:end));

