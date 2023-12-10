function J2=cost_second_opt(C)

global GM timestep tsol1 xsol1 dist_moon_min C3_min time_to_moon DV2 DV3

% Time span - limit of 100 days from apogee to reach lunar orbit
tspan2 = tsol1(end):timestep:tsol1(end)+100*24*3600;
% Apogee state vector after mid-course correction maneuver
CI0([1:3,13:15])   = xsol1(end,[1:3,13:15]) - xsol1(end,[10:12,22:24]) + ...
                     [0,0,0,C(1),C(2),C(3)];
CI0([4:6,16:18])   = xsol1(end,[4:6,16:18]) - xsol1(end,[10:12,22:24]);
CI0([7:9,19:21])   = xsol1(end,[7:9,19:21]) - xsol1(end,[10:12,22:24]);
CI0([10:12,22:24]) = [0,0,0,0,0,0];
% Integration
fbpfun = fourbp(GM);
[tsol2, xsol2] = ode89(fbpfun,tspan2,CI0);

% Minimum distance from Moon 
dist_moon = zeros(1,length(tsol2));
for j=1:length(tsol2)
    dist_moon(j) = norm(xsol2(j,1:3) - xsol2(j,7:9));
end
[dist_moon_min,ind_dist_min] = min(dist_moon);
% Time to arrival in Moon orbit
time_to_moon = tsol2(ind_dist_min);

% Minimum C3
C3 = zeros(1,length(tsol2));
for j=1:length(tsol2)
    C3(j) = ...
        norm(xsol2(j,13:15)-xsol2(j,19:21))^2 - 2*GM(3)/norm(xsol2(j,1:3)-xsol2(j,7:9));
end
[C3_min,~] = min(C3);

%%% COST FUNCTION %%%
% Correction maneuver Delta V
DV2  = norm(C(1:3));
% Lunar insertion Delta V
DV3  = norm(xsol2(ind_dist_min,13:15)-xsol2(ind_dist_min,19:21)) - sqrt(GM(3)/dist_moon_min);
% Cost function
J2 = abs(DV2)+abs(DV3);


