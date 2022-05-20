%% Call this script to postprocess all output .h5 files from one RunJobs sh.
clc,clear;
close all;

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
drs = [3, 6]; 
load = [5];
Vw = [2];
fw = [0.2];
theta1 = [0.0064];
theta2 = [-9];
A = [0.008];
AmB = [0.005];
NULoad = [0];
stress_dist = [1, 2, 3];
mesh = [2];
z_location = [0.000, 0.003, 0.004, 0.005];

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
Distance_To_Surface = 0.005;
for i = 1:1:size(filenames, 1)
    disp(filenames(i));
    WriteFiguresTet_function(filenames(i), Distance_To_Surface);
    close all;
    
    for iStress_dist = 1:1:size(stress_dist, 2)
        FakeDICDisp_function(filenames(i), stress_dist(iStress_dist));
        close all;
        
        RealStressAtLocation_function(filenames(i), stress_dist(iStress_dist), 0.005);
        close all;
        
    end
    
    for iZ_location = 1:1:size(z_location, 2)
        SlipRateAtDistInTheFault_function(filenames(i), z_location(iZ_location));
        close all;
    end
    
    % close all;
end