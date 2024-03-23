% This script calls fV_vs_slip_time function
clc,clear;
close all;
totalprefix = 'OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0';
% totalprefix = 'VaryingBOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005_0.006bD2000000.0Load5_Vw2_fw0.1_theta0.006_-9_NULoad2dir0_duration200_0';
totalprefix = 'OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.006_-9_NULoad2dir0_duration200_0';
totalprefix = 'OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0';
% totalprefix = "Hom0.0015_0.47_SW_0.1_2_36_40_W8_4ExcCorr2NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load8_Vw2_fw0.1_theta0.2208_-9_NULoad2dir0_duration200_50";
z_locations = [0, 5];
target_x = [10, 15, 20];

for iZ = 1:1:size(z_locations, 2)
    fV_vs_slip_time_function(totalprefix, target_x, z_locations(iZ), 0);
    close all;
end