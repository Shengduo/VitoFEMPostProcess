clc,clear;
close all;
videoprefix = '1WithWallDRS1.5_1.5ModA0.008Load5_Vw2_fw0.1_theta0.036_8000000000000000.0_NULoad2dir0';
dataFile = strcat('../dumpFiles/', videoprefix, '-domain.h5');

% Calculate normal direction
faultST = [-0.100000000000000, -0.0554300000000000];
faultND = [0.100000000000000, 0.0554300000000000];
normal_dir = (faultND - faultST);
normal_dir = normal_dir ./ norm(normal_dir, 2);
tangential_dir = [normal_dir, 0];
normal_dir = [-normal_dir(2), normal_dir(1), 0.];

% Wire and VS regions
WirePos1 = [-0.025657; -0.014222; 0];
VSstart = [0.006354, 0.003522, 0]';
VSend = [0.063204, 0.035034, 0]';
SmallVwStart = VSstart + tangential_dir' * 15 * 1e-3;
SmallVwEnd = SmallVwStart + tangential_dir' * 2 * 1e-3;
wireLabel = strcat('X = 0 [mm], wire position');
VSLabel = 'VS Region [mm]';

%% Process the domain file for displacements
h5disp(dataFile);
time = h5read(dataFile, '/time');
time = reshape(time, [1, size(time, 3)]);
displacement = h5read(dataFile, '/vertex_fields/displacement');
velocity = h5read(dataFile, '/vertex_fields/velocity');
Connection = h5read(dataFile, '/topology/cells');
Connection = Connection + 1;
nodalXYZ = h5read(dataFile, '/geometry/vertices');

% Select only the front surface
% Only extract the surface nodes
surface = 0.005;
NID = (abs(nodalXYZ(3, :) - surface) <= 1e-7);
EID = logical(zeros(1, size(Connection, 2)));
for i = 1:1:size(EID, 2)
    EID(i) = (sum(NID(Connection(:, i)')) == 3);
    if (EID(i) == 1)
        Connection(1:3, i) = Connection(NID(Connection(:, i)) == 1, i);
    end
end

% Update nodes and collection
NXYZs = nodalXYZ(1:2, NID);
nodeID = 0;
NIDD = zeros(1, size(NID, 2));
for i = 1:1:size(NID, 2)
    if NID(i) == 1
        nodeID = nodeID + 1;
        NIDD(i) = nodeID;
    end
end

Eles = Connection(1:3, EID);
for i = 1:1:size(Eles, 2)
    Eles(:, i) = NIDD(Eles(:, i))';
end

displacement = displacement(1:2, NID, :);
velocity = velocity(1:2, NID, :);

nOfTimeSteps = size(time, 2);

% Apply rotation to displacement and nodal XYZ to align with the 
% fault coordinate system
Q = [tangential_dir(1:2); normal_dir(1:2)];

% Rotate the velocity field
for i = 1:1:size(velocity, 3)
    velocity(:, :, i) = Q * velocity(:, :, i);
    displacement(:, :, i) = Q * displacement(:, :, i);
end

% Rotate nodal XYZ
XYZ = Q * NXYZs;

%% Get the interpolation of slip at these locations
% Get the start and end location of  VS region
VS_start_x = norm(VSstart - WirePos1, 2);
VS_end_x = norm(VSend - WirePos1, 2);

% Camera specific values
Pix_10mm = 81;  %Number of pixels in a grid of 10 mm ;  144 pixels in 2 mm -> 72 pixels in 1 mm
pxsize=10/(Pix_10mm)*1e-3; % pixel size in m/pixel  (Pix_10mm pixels  in 10 mm)
stepsize = 1;

% First specific the x and y's
x_up = VS_start_x : stepsize * pxsize : VS_end_x;
x_low = x_up;
y_up = pxsize : stepsize * pxsize : 20e-3;     % Assume the window height is 20 mm
y_low = - y_up;

% Get the values of displacements at [x, y]s

