close all
clear

tic

%%% GLOBAL VARIABLES AND KERNEL INITIALIZATION %%%
global GM  DEP_oe timestep rad_arr dist_earth_max C3_min dist_moon_min
global time2max_dist tsol1 xsol1 time_to_moon

% Load kernels
cspice_furnsh( 'metakr.tm' );

%--------------------------------------------------------------------------
%%% INITIAL GUESS & USER-DEFINED PARAMETERS %%%
% Launch date in UTC
timstr = '2024 JAN 21 05:00:00';  
% Spacecraft mass
mass_sc = 1000; %[kg]

% Departure ecliptic orbital elements (as used in cspice_conics and oscelt)
DEP_oe(1) = 6371 + 250; % radius of orbit in km
DEP_oe(2) = 0.9882036099;  % eccentricity
DEP_oe(3) = 0.3494;     % inclination
DEP_oe(4) = 0.1794;     % longitude of ascending node
DEP_oe(5) = 1.7197;     % argument of pericenter
% Mean anomaly - DO NOT CHANGE
DEP_oe(6) = 0;          

% Desired arrival radius (inclination to be added at a later time)
rad_arr = 1737 + 100;   % radius of arrival moon orbit in km

% Integration timestep. A step of 60 seconds is recommended for acceptable 
% computation times and accuracy.
timestep = 60; 

%--------------------------------------------------------------------------
%%% PARAMETERS %%%
% Launch date in seconds past J2000
t0 = cspice_str2et(timstr);  

%%% Gravitational parameters %%%
GM    = zeros(1,4);
% Spacecraft parameter
GM(1) = mass_sc * 6.67139e-20; 
% Celestial bodies parameters
GM(2) = cspice_bodvrd('SUN', 'GM', 1);
GM(3) = cspice_bodvrd('MOON', 'GM', 1);
GM(4) = cspice_bodvrd('EARTH', 'GM', 1);

%--------------------------------------------------------------------------
%%% FIRST OPTIMIZATION %%%
% Initial guess vector 
S0(1) = DEP_oe(2);
S0(2) = DEP_oe(3);
S0(3) = DEP_oe(4);
S0(4) = DEP_oe(5);
S0(5) = t0;
% Lower and upper bounds
lb = [0.975, 0,    0,    0,    t0 - 10*24*3600];
ub = [1,     pi/2, 2*pi, 2*pi, t0 + 10*24*3600];
% Options
options = optimoptions("fmincon",...
    "Algorithm","active-set",...
    "EnableFeasibilityMode",true,...
    "ConstraintTolerance",1e-05,...
    'Display','iter',...
    "MaxFunctionEvaluations",500);

% Plot initial guess and manually check quality. Execution can be stopped
% here if initial conditions are not satisfying.
%
% [Substitute at a later time with function that automatically checks 
% initial Earth-Sun-Moon geometry and add initial parameter optimization.]
CheckInitialGuess(S0)

% Continue or stop execution prompt
kbhit = input('\nType y to continue, anything else to stop: ','s');
if kbhit~='y'
    close all
    cspice_kclear
    return
else
    kbhit = [];
    close all
end

% First optimization execution
[S0_opt, DV13_opt] = ...
    fmincon(@cost_first_opt,S0,[],[],[],[],lb,ub,@constraints_first_opt,options);

%%% Recomputation of trajectory to apogee %%%
tspan1 = S0_opt(5):timestep:S0_opt(5)+time2max_dist;
% Transformation of orbital elements to cartesian state vector in ECI
DEP_cart = cspice_conics([DEP_oe,S0_opt(5),GM(4)]',S0_opt(5));
% Departure state vector
SI0([1:3,13:15])      = DEP_cart(1:6);
[SI0([4:6,16:18]), ~] = cspice_spkezr('sun',S0_opt(5),'j2000','NONE','earth');
[SI0([7:9,19:21]), ~] = cspice_spkezr('moon',S0_opt(5),'j2000','NONE','earth');
SI0([10:12,22:24])    = [0,0,0,0,0,0];
% Integration
fbpfun = fourbp(GM);
[tsol1, xsol1] = ode89(fbpfun,tspan1,SI0);


%--------------------------------------------------------------------------
%%% SECOND OPTIMIZATION %%%
% Initial guess - correction maneuver velocity components
C0  = [0,0,0];
% Lower and upper bounds
lb2 = [-1, -1, -1];
ub2 = [ 1,  1,  1];
% Second optimization execution
[C0_opt, DV23_opt] = ...
    fmincon(@cost_second_opt,C0,[],[],[],[],lb2,ub2,@constraints_second_opt,options);

%%% Recomputation of trajectory from apogee to Moon %%%
tspan2 = tsol1(end):timestep:time_to_moon;
% Apogee state vector after mid-course correction maneuver
CI0([1:3,13:15])   = xsol1(end,[1:3,13:15]) - xsol1(end,[10:12,22:24]) + ...
                     [0,0,0,C0_opt(1),C0_opt(2),C0_opt(3)];
CI0([4:6,16:18])   = xsol1(end,[4:6,16:18]) - xsol1(end,[10:12,22:24]);
CI0([7:9,19:21])   = xsol1(end,[7:9,19:21]) - xsol1(end,[10:12,22:24]);
CI0([10:12,22:24]) = [0,0,0,0,0,0];
% Integration
fbpfun = fourbp(GM);
[tsol2, xsol2] = ode89(fbpfun,tspan2,CI0);


%--------------------------------------------------------------------------
%%% PLOT %%%
Xsc1 = xsol1(:,1:3) - xsol1(:,10:12);
Xsc2 = xsol2(:,1:3) - xsol2(:,10:12);
Xm1  = xsol1(:,7:9) - xsol1(:,10:12);
Xm2  = xsol2(:,7:9) - xsol2(:,10:12);
plot3(Xsc1(:,1),Xsc1(:,2),Xsc1(:,3),'b')
hold on
plot3(Xsc2(:,1),Xsc2(:,2),Xsc2(:,3),'b')
plot3(Xm1(:,1),Xm1(:,2),Xm1(:,3),'r')
plot3(Xm2(:,1),Xm2(:,2),Xm2(:,3),'r')
% Draw Earth
[x,y,z] = sphere;
x = x*6371;
y = y*6371;
z = z*6371;
surf(x,y,z, 'FaceColor', '#4DBEEE')
% Draw Moon at lunar insertion
[x2, y2, z2] = sphere;
x2 = x2 * 1737 + Xm2(end,1);
y2 = y2 * 1737 + Xm2(end,2);
z2 = z2 * 1737 + Xm2(end,3);
surf(x2,y2,z2, 'FaceColor', '#808080')
% Options
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
grid on
axis equal
hold off

% Print results
elapsedTime = toc/60;
ToF = (tsol2(end) - S0_opt(5))/(24*3600); 
DepDate = cspice_et2utc(S0_opt(5),'C',0);
fprintf('\nDeparture date: %s\n', DepDate)
fprintf('Total Delta V: %.10f km/s\n', DV1+DV2+DV3)
fprintf('Departure from parking orbit - Delta V 1: %.5f km/s\n', DV1)
fprintf('Correction maneuver at apogee - Delta V 2: %.5f km/s\n', DV2)
fprintf('Lunar orbit insertion - Delta V 3: %.5f km/s\n', DV3)
fprintf('Time of Flight: %.2f days\n', ToF)
fprintf('Elapsed time: %.2f minutes\n', elapsedTime)

% Unload kernels
cspice_kclear
