%% Call this script to postprocess all output .h5 files from one RunJobs sh.
clc,clear;
close all;

%% First write figures for the default case
% filename = "1WithWallDRS1.5_1.5ModA0.008AmB0.005Load5_Vw2e+16_fw0.58_theta0.036_-11_NULoad2dir0";
filename = "1WithWallDRS1.5_1.5ModA0.008AmB0.005Load5_Vw2_fw0.1_theta0.036_-11_NULoad2dir0";
disp(filename);
% WriteFiguresTet_function(filenames(i), Distance_To_Surface);
% close all;
FakeDICDisp_function(filename);
close all;


%% Then all grid cases
% Parameters
drs = [1.5, 3]; 
load = [5];
Vw = [2];
fw = [0.1];
theta1 = [0.036];
theta2 = [-11, -10];
A = [0.008, 0.012, 0.016];
AmB = [0.005, 0.010];
NULoad = [0];

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
                                    fileNamePrefix = "1WithWallDRS1.5_" + string(drs(idrs)) + "ModA" + string(A(iA)) + "AmB" + string(AmB(iAmB)) + "Load" + string(load(iload)) + "_Vw" + string(Vw(iVw)) + "_fw" + string(fw(ifw)) + "_theta" + string(theta1(itheta1)) + "_" + string(theta2(itheta2)) + "_NULoad2dir" + string(NULoad(iNULoad));
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

%% Call the function script to postprocess all the cases in on RunJobs

Distance_To_Surface = 0.005;
for i = 1:1:size(filenames, 1)
    disp(filenames(i));
    % WriteFiguresTet_function(filenames(i), Distance_To_Surface);
    % close all;
    FakeDICDisp_function(filenames(i));
    close all;
end