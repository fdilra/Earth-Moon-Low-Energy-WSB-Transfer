function CheckInitialGuess(S)

global GM timestep tf DEP_oe

tspan = S(5):timestep:S(5)+150*24*3600;

% Transform departure orbital elements into cartesian state vector
DEP_cart = cspice_conics([DEP_oe, S(5), GM(4)]',S(5));
% Departure state vector of four bodies
SI0([1:3,13:15]) = DEP_cart(1:6);
[SI0([4:6,16:18]), ~] = cspice_spkezr('sun',S(5),'j2000','NONE','earth');
[SI0([7:9,19:21]), ~] = cspice_spkezr('moon',S(5),'j2000','NONE','earth');
SI0([10:12,22:24])   = [0,0,0,0,0,0];

% Integration
fbpfun = fourbp(GM);
[tsol, xsol] = ode89(fbpfun,tspan,SI0);

%--------------------------------------------------------------------------
% Max distance from Earth - Apogee
dist_earth = zeros(1,length(tsol));
for j=1:length(tsol)
    dist_earth(j) = norm(xsol(j,1:3) - xsol(j,10:12));
end
[~,ind_max_dist] = max(dist_earth);

% Minimum distance from Moon in second part of trajectory
dist_moon = zeros(1,length(tsol));
for j=1:length(tsol)
    dist_moon(j) = norm(xsol(j,1:3) - xsol(j,7:9));
end
dist_moon = dist_moon(ind_max_dist:end);
[dist_min_moon,~] = min(dist_moon);

% Characteristic energy in second part of trajectory
for j=1:length(tsol)
    C3(j) = ...
        norm(xsol(j,13:15)-xsol(j,19:21))^2 - 2*GM(3)/norm(xsol(j,1:3)-xsol(j,7:9));
end
C3 = C3(ind_max_dist:end);
[C3_min,~] = min(C3);

% Display minimum distance from Moon and minimum characteristic energy
fprintf('Minimum C3 = %.5f km^2/s^2\n', C3_min)
fprintf('Minimum distance from Moon: %.2f km\n', dist_min_moon)

%--------------------------------------------------------------------------
% Coordinates wrt ECI %
xSC  = xsol(:,1:3) - xsol(:,10:12);
xM   = xsol(:,7:9) - xsol(:,10:12);
% Sun position is scaled to be easily compared with the Moon and spacecraft
% positions at departure
xS   = (1/100)*(xsol(:,4:6) - xsol(:,10:12));

% Spacecraft trajectory
plot3(xSC(1:15:end,1),xSC(1:15:end,2),xSC(1:15:end,3))
hold on

% Draw Earth
[x,y,z] = sphere;
x = x*6371;
y = y*6371;
z = z*6371;
surf(x,y,z, 'FaceColor', '#4DBEEE')

% Moon orbit
plot3(xM(1:15:end,1),xM(1:15:end,2),xM(1:15:end,3))
% Draw Moon at departure
[x2, y2, z2] = sphere;
x2 = x2 * 1737 + xM(1,1);
y2 = y2 * 1737 + xM(1,2);
z2 = z2 * 1737 + xM(1,3);
surf(x2,y2,z2, 'FaceColor', '#FF0000')
% Draw Moon at closest approach
[x2, y2, z2] = sphere;
x2 = x2 * 1737 + xM(end,1);
y2 = y2 * 1737 + xM(end,2);
z2 = z2 * 1737 + xM(end,3);
surf(x2,y2,z2, 'FaceColor', '#808080')

% Draw Sun at departure
[x3, y3, z3] = sphere;
x3 = x3 * 10000 + xS(1,1);
y3 = y3 * 10000 + xS(1,2);
z3 = z3 * 10000 + xS(1,3);
surf(x3,y3,z3, 'FaceColor', 'y')

% Options
grid on
axis equal
hold off


