function [c,ceq]=constraints_second_opt(C)

global dist_earth_max C3_min rad_arr dist_moon_min time2max_dist

c(1) = C3_min + 0.0001;
c(2) = time2max_dist/(24*3600) - 50;
c = c';

ceq = dist_moon_min - rad_arr;
