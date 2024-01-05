function [s1,r1,ee1] = massBalance_Evap( s, u, q, moy ,discr_s, sys_param)
% Output:
%      s1 - final storage. 
%      r1 - release over the 24 hours. 
%
% Input:
%       s - initial storage. 
%       u - release decision 
%       q - inflow 
%       dt - day of the year
%
% See also MIN_RELEASE, MAX_RELEASE

HH = sys_param.integration_substep;
delta = sys_param.simulation.delta/HH; % Y
s_ = nan(HH+1,1);
r_ = nan(HH+1,1);
evap = interp1(discr_s, sys_param.simulation.ev(:,moy),s,'pchip',max(sys_param.simulation.ev(:,moy))); % evap in MCM/Y: lookup expected evap based on initial storage
%evap = interp1(discr_s, sys_param.simulation.ev(:,moy),s,'linear','extrap');

s_(1) = s;
s_max = discr_s(end); % maximum effective reservoir storage

for i=1:HH
    s_eff = max(s_(i) + delta*q - delta*evap,0); % effective storage
    [r_min, r_max] = min_max_relG3(s_eff, moy, s_max, sys_param) ;
    r_(i+1) = min( r_max , max( r_min , u ) ); % MCM/Y
    s_(i+1) = s_eff - delta*r_(i+1); % MCM
end

s1 = s_(HH);
r1 = mean(r_(2:end));
ee1 = evap;

