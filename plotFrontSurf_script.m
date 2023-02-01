%% Call this script to postprocess all output .h5 files from one RunJobs sh.
setEnvironment;

%% First write figures for the default case
filenames = ["3OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.006_-9_NULoad2dir0_duration200_0", ...
             "3OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.006_-6_NULoad2dir0_duration200_0", ...
             "GougeThick2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.006_-6_NULoad2dir0_duration200_0"];


times = [20, 32, 40, 48, 59, 68, 80, 88, 92];
times = times(2);
filenames = filenames(2);
% Write FrontFiles
for i = 1 : 1 : length(filenames)
    disp(filenames(i));
    
    % Write frontfile figures
    plotFrontSurf_function(filenames(i), times);
    close all;
    
end