function J1 = cost_first_opt(S)

global GM timestep DEP_oe dist_earth_max C3_min dist_moon_min time2max_dist
global DV1

% Time span
tspan = S(5):timestep:S(5)+150*24*3600;
% Transformation of orbital elements to cartesian state vector in ECI
DEP_cart = cspice_conics([DEP_oe,S(5),GM(4)]',S(5));
% Departure state vector
SI0([1:3,13:15])      = DEP_cart(1:6);
[SI0([4:6,16:18]), ~] = cspice_spkezr('sun',S(5),'j2000','NONE','earth');
[SI0([7:9,19:21]), ~] = cspice_spkezr('moon',S(5),'j2000','NONE','earth');
SI0([10:12,22:24])    = [0,0,0,0,0,0];
% Integration
fbpfun = fourbp(GM);
[tsol, xsol] = ode89(fbpfun,tspan,SI0);

% Maximum distance from Earth/Apogee determination
dist_earth = zeros(1,length(tsol));
for j=1:length(tsol)
    dist_earth(j) = norm(xsol(j,1:3) - xsol(j,10:12));
end
[dist_earth_max,ind_dist_max] = max(dist_earth);
time2max_dist = tsol(ind_dist_max)-S(5);

% Minimum distance from Moon in second arc
dist_moon = zeros(1,length(tsol));
for j=1:length(tsol)
    dist_moon(j) = norm(xsol(j,1:3) - xsol(j,7:9));
end
dist_moon = dist_moon(ind_dist_max:end);
[dist_moon_min,ind_dist_min] = min(dist_moon);
% Error prevention in case apogee is reached at time limit, in which case
% ind_dist_max+ind_dist_min exceeds array bounds in DV3.
if ind_dist_max == length(tsol)
    ind_dist_min = 0;
end

% Minimum C3 wrt Moon in second arc
C3 = zeros(1,length(tsol));
for j=1:length(tsol)
    C3(j) = ...
        norm(xsol(j,13:15)-xsol(j,19:21))^2 - 2*GM(3)/norm(xsol(j,1:3)-xsol(j,7:9));
end
C3 = C3(ind_dist_max:end);
[C3_min,~] = min(C3);

%%% COST FUNCTION %%%
% Departure Delta V
DV1  = abs(norm(xsol(1,13:15)) - sqrt(GM(4)/DEP_oe(1)));
% Lunar circular orbit insertion Delta V
DV3  = norm(xsol(ind_dist_max+ind_dist_min,13:15) - ...
    xsol(ind_dist_max+ind_dist_min,19:21)) - sqrt(GM(3)/dist_moon_min);
% Cost function 
J1 = abs(DV1) + abs(DV3);

