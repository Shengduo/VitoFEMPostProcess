% Read results from hdf5 files.
clc,clear;
close all;
videoprefix = '1WithWallDRS1.5_1.5ModA0.008AmB0.005Load5_Vw2e+16_fw0.58_theta0.036_-11_NULoad2dir0';
faultFileName = strcat('../faultFiles/', videoprefix, '-fault.h5');
frontsurfFile = strcat('../frontsurfFiles/', videoprefix, '-frontsurf.h5');

h5disp(faultFileName);
h5disp(frontsurfFile);
fontsize = 25;

% Read time
time = h5read(faultFileName, '/time');
time = reshape(time, [1, size(time, 3)]);
nOfTimeSteps = size(time, 2);

% Read node geometry
nodalXYZ = h5read(faultFileName, '/geometry/vertices');

% Select the surface nodes
I1 = (abs(nodalXYZ(3, :) - 0.005) <= 1e-7);
nOfNodes = size(nodalXYZ, 2);

% Read nodal slip, slip rate and traction
Slip = h5read(faultFileName, '/vertex_fields/slip');
SlipRate = h5read(faultFileName, '/vertex_fields/slip_rate');
Traction = h5read(faultFileName, '/vertex_fields/traction');
% Connection = h5read(faultFileName, '/topology/cells');

% Read nodal slip
% Input the first wire position
% WirePos1 = [-0.024138; -0.013380; 0];
wireLabel = strcat('X = 0 [mm], wire position');
% VSstart = [0.020468, 0.011345, 0]';
% VSend = [0.079942, 0.044312, 0]';
VSLabel = 'VS Region [mm]';

% Input the first wire position
WirePos1 = [-0.025657; -0.014222; 0];
% VSstart = [0.020468, 0.011345, 0]';
% VSend = [0.079942, 0.044312, 0]';

VSstart = [0.006354, 0.003522, 0]';
VSend = [0.063204, 0.035034, 0]';
VSregion = 1e3 * [norm(VSstart - WirePos1, 2), norm(VSend - WirePos1, 2)];

nOf2DNodes = sum(I1);

FaultX = zeros(1, nOf2DNodes);
surfaceNodesXYZ = nodalXYZ(:, I1);
[~, I] = sort(surfaceNodesXYZ(1,:));

% Fault Range
xrange = [0, 150];
yrange = [0, 10];

% Framerate for videos
framerate = 8;

for i = 1:1:nOf2DNodes
    FaultX(i) = norm(surfaceNodesXYZ(1:2, i) - WirePos1(1:2), 2) * sign(surfaceNodesXYZ(1, i) - WirePos1(1));
end
FaultX = FaultX(I);
surfaceNodesXYZ = surfaceNodesXYZ(1:2, I);

% 2D node coordinates
nodalXYZ2D = zeros(2, nOfNodes);
for i = 1:1:nOfNodes
    nodalXYZ2D(1, i) = sign(nodalXYZ(1, i) - WirePos1(1)) ...
        * norm(nodalXYZ(1:2, i) - WirePos1(1:2), 2);
    nodalXYZ2D(2, i) = nodalXYZ(3, i);
end

% Magnitude of slip rate
surfaceSlipRateMag = zeros(nOf2DNodes, nOfTimeSteps);
surfaceSlipRate = zeros(3, nOf2DNodes, nOfTimeSteps);

% Magnitude of slip
surfaceSlipMag = zeros(nOf2DNodes, nOfTimeSteps);
surfaceSlip = zeros(3, nOf2DNodes, nOfTimeSteps);

% Shear and Normal stress
surfaceShearStress = zeros(nOf2DNodes, nOfTimeSteps);
surfaceNormalStress = zeros(nOf2DNodes, nOfTimeSteps);

for t = 1:1:nOfTimeSteps
    surfaceSlipRate(:, :, t) = SlipRate(:, I1, t);
    surfaceSlip(:, :, t) = Slip(:, I1, t);
    surfaceShearStress(:, t) = Traction(1, I1, t);
    surfaceNormalStress(:, t) = Traction(3, I1, t);
    
    for i = 1:1:nOf2DNodes
        surfaceSlipRateMag(i, t) = norm(surfaceSlipRate(:, i, t));
        surfaceSlipMag(i, t) = norm(surfaceSlip(:, i, t));
    end
    surfaceSlipRateMag(:, t) = surfaceSlipRateMag(I, t);
    surfaceSlipMag(:, t) = surfaceSlipMag(I, t);
    surfaceShearStress(:, t) = surfaceShearStress(I, t);
    surfaceNormalStress(:, t) = surfaceNormalStress(I, t);
end

%% Plot V-Slip history at X [mm]
target_x = 60;

% Find the target index
[~, Ind] = min(abs(FaultX - target_x / 1e3));
targetCoords = surfaceNodesXYZ(:, Ind);

% Generate sliprate-slip plot
fig1 = figure(1);
pos = fig1.Position;
pos(3) = pos(3) * 3;
pos(4) = pos(4) * 2;
fig1.Position = pos;

subplot(2, 3, 4);
plot(surfaceSlipMag(Ind, 2:end) * 1e6, surfaceSlipRateMag(Ind, 2:end), 'linewidth', 2.0);
hold on; grid on;
scatter(surfaceSlipMag(Ind, 2:end) * 1e6, surfaceSlipRateMag(Ind, 2:end), 'filled');
xlabel('Slip [\mu m]');
ylabel('Slip rate [m/s]');
xlim(xrange);
ylim(yrange);
set(gca, 'FontSize', fontsize);

%% Plot friction-shearstress history at X [mm]
% target_x = 80;

% Find the target index
[val, Ind] = min(abs(FaultX - target_x / 1e3));

% Generate sliprate-slip plot
subplot(2, 3, 1);
plot(surfaceSlipMag(Ind, 2:end) * 1e6, - surfaceShearStress(Ind, 2:end) ./ surfaceNormalStress(Ind, 2:end), 'linewidth', 2.0);
hold on; grid on;
scatter(surfaceSlipMag(Ind, 2:end) * 1e6, - surfaceShearStress(Ind, 2:end) ./ surfaceNormalStress(Ind, 2:end), 'filled');
% xlabel('Slip [\mu m]');
ylabel('Friction coefficient');
xlim(xrange);
ylim([0, 1]);
probeLabel = strcat('X ={ }', num2str(target_x), '{ }[mm]');
title(probeLabel, 'Fontsize', fontsize);
set(gca, 'FontSize', fontsize);

%% Save a video of fault-normal velocity field
% Read time
time = h5read(frontsurfFile, '/time');
time = reshape(time, [1, size(time, 3)]);
nOfTimeSteps = size(time, 2);

% Read node geometry
nodalXYZ = h5read(frontsurfFile, '/geometry/vertices');
nodalXYZ = nodalXYZ(1:2, :);
nOfNodes = size(nodalXYZ, 2);

% Read nodal slip, slip rate and traction
velocity = h5read(frontsurfFile, '/vertex_fields/velocity');
velocity = velocity(1:2, :, :);
velocityMag = zeros(size(velocity, 2), size(velocity, 3));
for i = 1:1:size(velocity, 2)
    for j = 1:1:size(velocity, 3)
        velocityMag(i, j) = norm(velocity(:, i, j), 2);
    end
end
Connection = h5read(frontsurfFile, '/topology/cells');

% Rotate into fault coordinate system
alpha = 29 / 180 * pi;
Q = [cos(alpha), sin(alpha); -sin(alpha), cos(alpha)];


%% -------------------------------------------------------------------------
% Save a video for fV vs slip
% Other constants
timeWindow = [10, 120];
videoflag = true;
videoXrange = [-110, 110];
videoYrange = [-110, 110];
crange = [0, 2];

if videoflag == true
    
    % Initialize names
    videoname = strcat(videoprefix,'X_', num2str(target_x), '_f_and_V_vs_slip.avi');
    
    % Initialize video
    myVideo = VideoWriter(strcat('../Videos/', videoname), 'Motion JPEG AVI');
    myVideo.FrameRate = framerate;
    myVideo.Quality = 100;
    open(myVideo);
    
    % Shoot the video
    for i = 1:1:size(time, 2)
        if time(i) * 1e6 < timeWindow(1)
            continue;
        end
        if time(i) * 1e6 > timeWindow(2)
            break;
        end
        if exist('lastLine1', 'var')
            delete(lastLine1);
        end
        if exist('lastLine2', 'var')
            delete(lastLine2);
        end
        subplot(2, 3, 1);
        lastLine1 = xline(1e6 * surfaceSlipMag(Ind, i), 'linewidth', 2.0);
        hold on;
        subplot(2, 3, 4);
        lastLine2 = xline(1e6 * surfaceSlipMag(Ind, i), 'linewidth', 2.0);
        hold on;
        
        subplot(2, 3, [2, 3, 5, 6]);
        p = patch(1e3 * nodalXYZ(1, :), 1e3 * nodalXYZ(2, :), ...
            velocityMag(:, i), 'EdgeColor', 'none', 'FaceColor', 'interp');
        hold on;
        p.Faces = Connection' + 1;
        c = colorbar;
        caxis(crange);
        ylabel(c,'Particle Velocity Magnitude[m/s]','FontName','Avenir','FontSize',fontsize);
        axis equal;
        grid on;
        xlim(videoXrange);
        ylim(videoYrange);
        plot(1e3 * [surfaceNodesXYZ(1, 1), surfaceNodesXYZ(1, end)], ...
             1e3 * [surfaceNodesXYZ(2, 1), surfaceNodesXYZ(2, end)], ...
             'w', 'linewidth', 1.0);
        hold on;
        
        % Plot the VS region
        plot(1e3 * [VSstart(1), VSend(1)], ...
             1e3 * [VSstart(2), VSend(2)], ...
             'r', 'linewidth', 1.5);
        text(0, 60, VSLabel, 'color', 'red', 'Fontsize', fontsize);
        hold on;
        
        % Scatter the wire position
        scatter(1e3 * WirePos1(1), 1e3 * WirePos1(2), 80, 'filled', 'w');
        text(1e3 * WirePos1(1) + 5, 1e3 * WirePos1(2) - 2, wireLabel, 'color', 'w', 'Fontsize', fontsize);
        
        scatter(1e3 * targetCoords(1), 1e3 * targetCoords(2), 80, 'filled', 'y');
        text(1e3 * targetCoords(1) + 5, 1e3 * targetCoords(2) - 2, probeLabel, 'color', 'y', 'Fontsize', fontsize);
        hold off;
        % yline(0, 'LineWidth', 2.0, 'color', 'w');
        xlabel('Width [mm]');
        ylabel('Height [mm]');
        title(strcat('Time = ', num2str(1e6 * time(i) - 10, '%2f'), ' [\mu s]'));
        set(gca, 'FontSize', fontsize);
        
        % Write the video
        frame = getframe(gcf);
        writeVideo(myVideo, frame);
    end
    close(myVideo);
end

%% -------------------------------------------------------------------------
% Save a video for fV vs time
videoflag = true;
%% Plot friction-shearstress history at X [mm]
% target_x = 80;
% Generate sliprate-slip plot

fig2 = figure(2);
pos = fig2.Position;
pos(3) = pos(3) * 3;
pos(4) = pos(4) * 2;
fig2.Position = pos;
xrange = [0, 110];
% Generate sliprate-slip plot
subplot(2, 3, 1);
plot(time * 1e6 - 10, - surfaceShearStress(Ind, 1:end) ./ surfaceNormalStress(Ind, 2), 'linewidth', 2.0);
hold on; grid on;
scatter(time * 1e6 - 10, - surfaceShearStress(Ind, 1:end) ./ surfaceNormalStress(Ind, 2), 'filled');
% xlabel('Slip [\mu m]');
ylabel('Friction coefficient');
xlim(xrange);
ylim([0, 1]);
probeLabel = strcat('X ={ }', num2str(target_x), '{ }[mm]');
title(probeLabel, 'Fontsize', fontsize);
set(gca, 'FontSize', fontsize);

subplot(2, 3, 4);
plot(time * 1e6 - 10, surfaceSlipRateMag(Ind, 1:end), 'linewidth', 2.0);
hold on; grid on;
scatter(time * 1e6 - 10, surfaceSlipRateMag(Ind, 1:end), 'filled');
xlabel('Time [\mu s]');
ylabel('Slip rate [m/s]');
xlim(xrange);
ylim(yrange);
set(gca, 'FontSize', fontsize);


% Other constants
timeWindow = [10, 150];
% videoflag = true;
videoXrange = [-110, 110];
videoYrange = [-110, 110];
crange = [0, 2];

if videoflag == true
    
    % Initialize names
    videoname = strcat(videoprefix,'X_', num2str(target_x), '_f_and_V_vs_time.avi');
    
    % Initialize video
    myVideo = VideoWriter(strcat('../Videos/', videoname), 'Motion JPEG AVI');
    myVideo.FrameRate = framerate;
    myVideo.Quality = 100;
    open(myVideo);
    
    % Shoot the video
    for i = 1:1:size(time, 2)
        if time(i) * 1e6 < timeWindow(1)
            continue;
        end
        if time(i) * 1e6 > timeWindow(2)
            break;
        end
        if exist('lastLine1', 'var')
            delete(lastLine1);
        end
        if exist('lastLine2', 'var')
            delete(lastLine2);
        end
        subplot(2, 3, 1);
        lastLine1 = xline(1e6 * time(i) - 10, 'linewidth', 2.0);
        hold on;
        subplot(2, 3, 4);
        lastLine2 = xline(1e6 * time(i) - 10, 'linewidth', 2.0);
        hold on;
        
        subplot(2, 3, [2, 3, 5, 6]);
        p = patch(1e3 * nodalXYZ(1, :), 1e3 * nodalXYZ(2, :), ...
            velocityMag(:, i), 'EdgeColor', 'none', 'FaceColor', 'interp');
        hold on;
        p.Faces = Connection' + 1;
        c = colorbar;
        caxis(crange);
        ylabel(c,'Particle Velocity Magnitude[m/s]','FontName','Avenir','FontSize',fontsize);
        axis equal;
        grid on;
        xlim(videoXrange);
        ylim(videoYrange);
        plot(1e3 * [surfaceNodesXYZ(1, 1), surfaceNodesXYZ(1, end)], ...
             1e3 * [surfaceNodesXYZ(2, 1), surfaceNodesXYZ(2, end)], ...
             'w', 'linewidth', 1.0);
        hold on;
        
        % Plot the VS region
        plot(1e3 * [VSstart(1), VSend(1)], ...
             1e3 * [VSstart(2), VSend(2)], ...
             'r', 'linewidth', 1.5);
        text(0, 60, VSLabel, 'color', 'red', 'Fontsize', fontsize);
        hold on;
        
        % Scatter the wire position
        scatter(1e3 * WirePos1(1), 1e3 * WirePos1(2), 80, 'filled', 'w');
        text(1e3 * WirePos1(1) + 5, 1e3 * WirePos1(2) - 2, wireLabel, 'color', 'w', 'Fontsize', fontsize);
        
        scatter(1e3 * targetCoords(1), 1e3 * targetCoords(2), 80, 'filled', 'y');
        text(1e3 * targetCoords(1) + 5, 1e3 * targetCoords(2) - 2, probeLabel, 'color', 'y', 'Fontsize', fontsize);
        hold off;
        % yline(0, 'LineWidth', 2.0, 'color', 'w');
        xlabel('Width [mm]');
        ylabel('Height [mm]');
        title(strcat('Time = ', num2str(1e6 * time(i) - 10, '%2f'), ' [\mu s]'));
        set(gca, 'FontSize', fontsize);
        
        % Write the video
        frame = getframe(gcf);
        writeVideo(myVideo, frame);
    end
    close(myVideo);
end


