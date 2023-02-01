% This script calls fV_vs_slip_time function
clc,clear;
close all;
setEnvironment;
% totalprefix = 'OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0';
% totalprefix = 'VaryingBOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005_2.0bD10Load5_Vw2_fw0.1_theta0.006_-9_NULoad2dir0_duration200_0';
% totalprefix = 'VaryingBOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Lw3Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0';
% totalprefix = 'OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.006_-9_NULoad2dir0_duration200_0';
% totalprefix = 'OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.0045_-9_NULoad2dir1_duration200_0';
% totalprefix = 'VaryingBOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005_3.0bD30Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0';
% totalprefix = 'VaryingBOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005_2.0bD7Load5_Vw2_fw0.1_theta0.0045_-9_NULoad2dir1_duration200_0';
totalprefix = 'OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.006_-6_NULoad2dir0_duration200_0';
% totalprefix = 'VwSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_10Si11_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0';
z_locations = [0, 5];
target_x = [5, 15, 20];

for iZ = 1:1:size(z_locations, 2)
    fV_vs_slip_time_function_nomesh(totalprefix, target_x, z_locations(iZ), 0, 1);
    close all;
end