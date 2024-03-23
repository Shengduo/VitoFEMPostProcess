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
% totalprefix = 'OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.006_-6_NULoad2dir0_duration200_0';
% totalprefix = 'W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load8_Vw2_fw0.1_theta0.3399_-8_NULoad2dir0_duration200_8';
% totalprefix = 'VwSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_10Si11_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0';
% totalprefix = 'SW_0.1_0.01_W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load8_Vw2_fw0.1_theta0.3399_-8_NULoad2dir0_duration200_8';
% totalprefix = "Hom0.0015_0.47_SW_0.1_2_36_40_W8_4ExcCorr2NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load8_Vw2_fw0.1_theta0.2208_-9_NULoad2dir0_duration200_50";
% totalprefix = "Blk10_8_W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.012AmB0.01Load8_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8";
% totalprefix = "BlkW8_SW_0.01_2_8_8_36_40_7e-7_W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.012AmB0.01Load8_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_13";
% totalprefixs = ["W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8", ...
%                 "W8_4ExcCorr2NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8", ...
%                 "W8_4ExcCorr3NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8"];

totalprefixs = ["W8_4ExcCorr2NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8"];

% z_locations = [0, 5];
% target_x = [5, 15, 20];
% totalprefixs = ["SW_HomaGouge_W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.003AmB-0.005Load5_Vw0.2_fw0.33_theta0.1434_0.1434_NULoad2dir0_duration200_8"];
% totalprefixs = ["SlipW8_4ExcCorr2NPDirWithWallDRS0.8ModA0.016AmB0.014Load5_Vw2_fw0.1_theta-7_-7_NULoad2dir0_duration200_8"];
z_locations = [0];
target_x = [-5, 5];

for iZ = 1:1:size(z_locations, 2)
    fV_vs_slip_time_function_nomesh_allLabel(totalprefixs, target_x, z_locations(iZ), 0, 2);
    close all;
end