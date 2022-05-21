% This script calls fV_vs_slip_time function
clc,clear;
close all;
totalprefix = 'WithWallDRS1.5_1.5ModA0.008AmB0.005Load5_Vw2_fw0.1_theta0.036_-11_NULoad2dir0';
z_locations = [5, 0];
target_x = [37, 42, 57];

for iZ = 1:1:size(z_locations, 2)
    fV_vs_slip_time_function(totalprefix, target_x, z_locations(iZ));
    close all;
end