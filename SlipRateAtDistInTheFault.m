clc,clear;
% This script process the faultfile and get the slip at certain location

load('BRColorScale.mat');
VitoColorFlag = true;
tractionOffsetFlag = false;

videoprefix = '1WithWallDRS1.5_1.5ModA0.008AmB0.005Load5_Vw2_fw0.1_theta0.036_-11_NULoad2dir0';
faultFileName = strcat('../faultFiles/', videoprefix, '-fault.h5');
h5disp(faultFileName);
fontsize = 25;

% Read time
time = h5read(faultFileName, '/time');
time = reshape(time, [1, size(time, 3)]);
time = time - 10e-6;
nOfTimeSteps = size(time, 2);

% Several wave speed to show on the aplot
cp = 2662.4;
cs = 1279;
nu = 0.35;
cr = (0.874 + 0.196 * nu - 0.043 * nu^2 - 0.055 * nu^3) * cs;
cX = [60, 80];
crY = [10, (cX(2) - cX(1)) * 1e3 / cr + 10];
csY = [10, (cX(2) - cX(1)) * 1e3 / cs + 10];
cpY = [10, (cX(2) - cX(1)) * 1e3 / cp + 10];

% Read node geometry
nodalXYZ = h5read(faultFileName, '/geometry/vertices');