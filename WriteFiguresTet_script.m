%% Call this script to postprocess all output .h5 files from one RunJobs sh.
setEnvironment;

%% First write figures for the default case
% filename = "1WithWallDRS1.5_1.5ModA0.008AmB0.005Load5_Vw2e+16_fw0.58_theta0.036_-11_NULoad2dir0";
filename = "1DirWithWallDRS1.5_1.5ModA0.008AmB0.005Load5_Vw2_fw0.1_theta0.036_-11_NULoad2dir0_duration120";
disp(filename);
% WriteFiguresTet_function(filenames(i), Distance_To_Surface);
% close all;
% FakeDICDisp_function(filename);
close all;


%% Then all grid cases
% Parameters
drs = [1.5, 3, 6]; 
load = [5];
Vw = [2];
fw = [0.2];
theta1 = [0.0064];
theta2 = [-9];
A = [0.011];
AmB = [0.005];
NULoad = [5];
stress_dist = [1, 2, 3];
mesh = [1];
% z_location = [0.005, 0.003, 0.000];
% z_location = [0.000, 0.005];
z_location = [0.005];

% Filenames the things are stored
filenames = [];
for idrs = 1:1:size(drs, 2)
    for iload = 1:1:size(load, 2)
        for iVw = 1:1:size(Vw, 2)
            for ifw = 1:1:size(fw, 2)
                for itheta1 = 1:1:size(theta1, 2)
                    for itheta2 = 1:1:size(theta2, 2)
                        for iA = 1:1:size(A, 2)
                            for iAmB = 1:1:size(AmB, 2)
                                if (A(iA) - AmB(iAmB) <= 0) 
                                    continue;
                                end
                                for iNULoad = 1:1:size(NULoad, 2)
                                    for iMesh = 1:1:size(mesh, 2)
                                        % fileNamePrefix =  string(mesh(iMesh)) + "WithWallDRS1.5_" + string(drs(idrs)) + "ModA" + string(A(iA)) + "AmB" + string(AmB(iAmB)) + "Load" + string(load(iload)) + "_Vw" + string(Vw(iVw)) + "_fw" + string(fw(ifw)) + "_theta" + string(theta1(itheta1)) + "_" + string(theta2(itheta2)) + "_NULoad2dir" + string(NULoad(iNULoad));
                                        fileNamePrefix =  string(mesh(iMesh)) + "NPDirWithWallDRS1.5_" + string(drs(idrs)) + "ModA" + string(A(iA)) + "AmB" + string(AmB(iAmB)) + "Load" + string(load(iload)) + "_Vw" + string(Vw(iVw)) + "_fw" + string(fw(ifw)) + "_theta" + string(theta1(itheta1)) + "_" + string(theta2(itheta2)) + "_NULoad2dir" + string(NULoad(iNULoad)) + "_duration120";
                                        filenames = [filenames; fileNamePrefix];
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Call the function script to postprocess all the cases in on RunJobs

% 2 Distance_To_Surface = 0.004133;
% 1 Distance_To_Surface = 0.004135;
% 3 Distance_To_Surface = 0.004134;

% filenames = ["1NPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.00224_-9_NULoad2dir0_duration120",
%              "2NPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.0016_-9_NULoad2dir1_duration120_100",
%              "2NPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.07_-9_NULoad2dir2_duration120_100", 
%              "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration120_0", 
%              "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.0045_-9_NULoad2dir1_duration120_0", 
%              "2NPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.0016_-9_NULoad2dir3_duration120_100", 
%              "2NPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.0016_-9_NULoad2dir3_duration120_0"];
         
%% Longer versions of shits
filenames = ["2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.006_-9_NULoad2dir0_duration200_0", ...
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.006_-9_NULoad2dir0_duration200_100", ...
             "2NPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.0016_-9_NULoad2dir3_duration200_0", ...
             "2NPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.0016_-9_NULoad2dir3_duration200_100", ...
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.0045_-9_NULoad2dir1_duration200_0", ...
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.0045_-9_NULoad2dir1_duration200_100"];
         
filenames = [filenames, ...
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ...
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_100", ...
             "2NPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.0016_-9_NULoad2dir1_duration200_0", ...
             "2NPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.0016_-9_NULoad2dir1_duration200_100", ...
             "2NPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.00224_-9_NULoad2dir0_duration200_0", ...
             "2NPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.00224_-9_NULoad2dir0_duration200_100"];
         
% Baseline cases
filenames = [filenames, ...
            "2OPDirWithWallDRS1.5_1.5ModA0.011AmB-0.005Load5_Vw0.2_fw0.33_theta0.006_0.006_NULoad2dir0_duration200_0", ...
            "2NPDirWithWallDRS1.5_1.5ModA0.011AmB-0.005Load5_Vw1.1_fw0.27_theta0.00224_0.00224_NULoad2dir0_duration200_0", ...
            "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2e+16_fw0.1_theta0.006_-9_NULoad2dir0_duration200_0", ...
            "2NPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2e+16_fw0.1_theta0.00224_-9_NULoad2dir0_duration200_0"];
         
         
pre_times = [0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 0, 0, 0] * 1e-6;
filenames = filenames([1, 5, 7, 13, 15]);
pre_times = pre_times([1, 5, 7, 13, 15]);
% filenames = ["1VaryingBOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005_3.0bD30Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ...
%             "2VaryingBOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005_2.0bD7Load5_Vw2_fw0.1_theta0.0045_-9_NULoad2dir1_duration200_0"];
filenames = ["2VaryingBOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005_2.0bD10Load5_Vw2_fw0.1_theta0.006_-9_NULoad2dir0_duration200_0"];
filenames = ["2VaryingBOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Lw10Load5_Vw2_fw0.1_theta0.006_-9_NULoad2dir0_duration200_0"];
filenames = ["1VaryingSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_10Si7.5_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0"];
% filenames = ["1VaryingBOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Lw3Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0"];
% filenames = ["1VaryingSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_3Si7.5_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0"];
filenames = ["1VWSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_5Si11_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0"];
Distance_To_Surface = 0.005;
% pre_time = 0e-6;
filenames = ["2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.006_-9_NULoad2dir0_duration200_0", ...      % 1
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2e+16_fw0.1_theta0.006_-9_NULoad2dir0_duration200_0", ...  % 2
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.006_-6_NULoad2dir0_duration200_0", ...      % 3
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2e+16_fw0.1_theta0.006_-6_NULoad2dir0_duration200_0", ...  % 4
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.006_-7_NULoad2dir0_duration200_0", ...      % 5
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2e+16_fw0.1_theta0.006_-7_NULoad2dir0_duration200_0", ...  % 6
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.006_-8_NULoad2dir0_duration200_0", ...
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2e+16_fw0.1_theta0.006_-8_NULoad2dir0_duration200_0", ...
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB-0.005Load5_Vw0.2_fw0.33_theta0.006_0.006_NULoad2dir0_duration200_0", ...
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ...      % 10
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.005_-8_NULoad2dir3_duration200_0", ...      % 11
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.005_-7_NULoad2dir3_duration200_0", ...      % 12     
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.005_-6_NULoad2dir3_duration200_0", ...      % 13
             "1OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.005_-7.7_NULoad2dir3_duration200_0", ...    % 14
             "1OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.005_-7.5_NULoad2dir3_duration200_0", ...    % 15
             "1OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.005_-7.3_NULoad2dir3_duration200_0", ...
             "1OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.005_-8.5_NULoad2dir3_duration200_0", ...
             "1OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.005_-8.7_NULoad2dir3_duration200_0", ...
             "1OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.005_-8.3_NULoad2dir3_duration200_0", ...
             "1VarySiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_3Si7.5_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ...    % 20
             "1VwSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_10Si11_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ...
             "1VaryExpSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_3Si7.5_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ...
             "1VaryExpSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_3 100Si7.5_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ...
             "1VaryExpSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_3 20Si7.5_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ...
             "1VaryExpSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_3 40Si7.5_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ...   % 25
             "1VaryExpSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_3 60Si7.5_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ...
             "1VaryExpSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_3 80Si7.5_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ...
             "1VaryExpAllSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_3_10Si7.5_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ...
             "1VaryExpAllSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_3_20Si7.5_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ...
             "1VaryExpAllSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_3_40Si7.5_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ...   % 30
             "1VaryExpAllSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_3_50Si7.5_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ...   % 31
             "1VaryExpAllSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_3_80Si7.5_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ...   % 32
             "1VaryExpAllSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_3_200Si7.5_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ...  % 33
             "1VaryExpAllSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_3_200Si7.5_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ...  % 34
             "1VaryExpAllSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_3_300Si7.5_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ...  % 35
             "1VaryExpAllSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_3_500Si7.5_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ...  % 36
             "1VaryExpAllSiOPDirWithWallDRS1.5_1.5ModA0.011AmB0.005L10000_3_1000Si7.5_15Load5_Vw2_fw0.1_theta0.005_-9_NULoad2dir3_duration200_0", ... % 37
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.0045_-7_NULoad2dir1_duration200_0", ...    % 38
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.0045_-7.5_NULoad2dir1_duration200_0", ...    % 39
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.0045_-8_NULoad2dir1_duration200_0", ...    % 40
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.0045_-8.5_NULoad2dir1_duration200_0", ...    % 41
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.0045_-9_NULoad2dir1_duration200_0", ...    % 42
             "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.0045_-6_NULoad2dir1_duration200_0", ...
             "GougeThick2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.006_-6_NULoad2dir0_duration200_0"];
pre_times = zeros(1, length(filenames));
filenames = filenames([9]);

% New files 1208
filenames = ["Corr1OPDirWithWallDRS1.5_1.5ModA0.011AmB0.01Load5_Vw2_fw0.1_theta0.006_-6_NULoad2dir0_duration200_0", ...
             "Corr1OPDirWithWallDRS1.5_1.5ModA0.011AmB0.01Load5_Vw2_fw0.1_theta0.006_-5_NULoad2dir0_duration200_0", ...
             "Corr1OPDirWithWallDRS1.5_1.5ModA0.011AmB0.01Load5_Vw2_fw0.1_theta0.006_-4_NULoad2dir0_duration200_0"]; 

filenames = [filenames, ...
             "ExcCorr1OPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.006_-7_NULoad2dir0_duration200_0", ...
             "ExcCorr1OPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.006_-8_NULoad2dir0_duration200_0", ...
             "ExcCorr1OPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.006_-9_NULoad2dir0_duration200_0"];
filenames = ["HMonly1210", ...
             "HMonly1210_theta1.2261", ...
             "HMonly1210_theta0.78181", ...
             "HMonly1210_theta0.67412", ...
             "HMonly1210_theta0.43217", ... % TOO SLOW
             "HMonly1210_theta0.32131", ... % TOO FAST
             "HMonly1210_theta0.37264", ... % TOO FAST
             "HMonly1210_theta0.4013", ...  % TOO FAST
             "HMonly1210_theta0.41954", ... % A BIT FAST
             "HMonly1210_theta0.42581", ... % A Little bit fast
             "HMonly1210_theta0.42898"];

filenames = ["ExcCorr1OPDirWithWallDRS1.5_1.5ModA0.011AmB-0.005Load10_Vw0.2_fw0.33_theta1.4664_1.4664_NULoad2dir0_duration200_15", ...
             "ExcCorr1OPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load10_Vw2_fw0.1_theta1.4664_-7_NULoad2dir0_duration200_20", ...
             "W8_4ExcCorr1OPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta1.4664_-7_NULoad2dir0_duration200_20", ...
             "W8_4_Load8_HMonly1212_3_8_theta0.1434", ...
             "W8_4_Load6_6HMonly1212_3_8_theta0.1434", ...
             "W8_4_Load5_8HMonly1212_3_8_theta0.1434", ...
             "W8_4_Load5_6HMonly1212_3_8_theta0.1434", ...
             "W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8", ... % 8 ## Ones to be used ##
             "W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load8_Vw2_fw0.1_theta0.3399_-8_NULoad2dir0_duration200_8", ...
             "W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load8_Vw2_fw0.1_theta0.8061_-9_NULoad2dir0_duration200_8", ...
             "W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.011AmB0.009Load8_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8", ...
             "Visco1e40_W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8", ...
             "Elastic_W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8", ...
             "Viscos_W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8", ...
             "Elastic_lowCsExcCorr1NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8", ... % 15
             "Visco_MaxDt_ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8", ... % 16
             "ChangingE_W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load10_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8", ... % 17
             "L100_ChangingE_W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8", ... % 18
             "L50_200_ChangingE_W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8", ... % 19
             "L100_ChangingE_W8_4.5_1.5ModA0.016AmB0.014Load10_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8", ... % 20
             "L50_200_damping2.0_ChangingE_W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8", ... % 21
             "F0.33_0.2_L50_1000000_ChangingE_cp3500_W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8", ... % 22
             "FINITETHICKNESSGOUGE_W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8"]; % 23

z_location = [0.005, 0];
filenames = filenames([22]);
pre_times = [0.e-6];
% filenames = ["W16_4ExcCorr1OPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load10_Vw2_fw0.1_theta1.4664_-7_NULoad2dir0_duration200_5", ...
%              "W16_4ExcCorr1OPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load10_Vw2_fw0.1_theta1.4664_-7_NULoad2dir0_duration200_15", ...
%              "W16_4ExcCorr1OPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load10_Vw2_fw0.1_theta1.4664_-7_NULoad2dir0_duration200_20", ...
%              "W16_4ExcCorr1OPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load10_Vw2_fw0.1_theta1.4664_-7_NULoad2dir0_duration200_10"];
% % for i = 12:1:12
for i = 1:1:size(filenames, 2)
    disp(filenames(i));
    % WriteFiguresTet_function(filenames(i), Distance_To_Surface, pre_times(i));
    close all;
    
%     for iStress_dist = 1:1:size(stress_dist, 2)
%         FakeDICDisp_function(filenames(i), stress_dist(iStress_dist));
%         close all;
%         
%         RealStressAtLocation_function(filenames(i), stress_dist(iStress_dist), 0.005);
%         close all;
%         
%     end
    
    for iZ_location = 1:1:size(z_location, 2)
      SlipRateAtDistInTheFault_function(filenames(i), z_location(iZ_location), pre_times(i), 'png');
      close all;
    end
    
%     close all;
%     writeVideosTet_function(filenames(i), pre_times(i));
%     close all;
end