clc,clear;
close all;
videoprefix = 'DRS1_2.0LongLoad4.75_Vw1.5_fw0.16_theta0.08';
upperFile = strcat(videoprefix, '-upper_crust.h5');
lowerFile = strcat(videoprefix, '-lower_crust.h5');

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
VSend = [0.058832, 0.032610, 0]';
SmallVwStart = VSstart + tangential_dir' * 15 * 1e-3;
SmallVwEnd = SmallVwStart + tangential_dir' * 2 * 1e-3;
wireLabel = strcat('X = 0 [mm], wire position');
VSLabel = 'VS Region [mm]';

%% Process the lower .hdf5 file for stress
h5disp(lowerFile);
stress = h5read(lowerFile, '/cell_fields/stress');
Connection = h5read(lowerFile, '/topology/cells');
Connection = Connection + 1;
nodalXYZ = h5read(lowerFile, '/geometry/vertices');
time = h5read(lowerFile, '/time');
time = reshape(time, [1, size(time, 3)]);

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
NXYZs_lower = nodalXYZ(1:2, NID);
nodeID = 0;
NIDD = zeros(1, size(NID, 2));
for i = 1:1:size(NID, 2)
    if NID(i) == 1
        nodeID = nodeID + 1;
        NIDD(i) = nodeID;
    end
end

Eles_lower = Connection(1:3, EID);
for i = 1:1:size(Eles_lower, 2)
    Eles_lower(:, i) = NIDD(Eles_lower(:, i))';
end

stress = stress(:, EID, :);

nOfTimeSteps = size(stress, 3);
nOfCells = size(stress, 2);
shearStress_lower = zeros(nOfCells, nOfTimeSteps);
n = normal_dir;
t = [-n(2), n(1), 0.];

% Calculate the shear stress as <n \cdot \sigma \cdot t>
for j = 1:1:nOfTimeSteps
    for i = 1:1:nOfCells
        sigma = [stress(1, i, j), stress(2, i, j), stress(3, i, j); ...
                 stress(2, i, j), stress(4, i, j), stress(5, i, j); ...
                 stress(3, i, j), stress(5, i, j), stress(6, i, j)];
        shearStress_lower(i, j) = n * sigma * t';
    end
end

%% Process the upper .hdf5 file for stress
h5disp(upperFile);
stress = h5read(upperFile, '/cell_fields/stress');
Connection = h5read(upperFile, '/topology/cells');
Connection = Connection + 1;
nodalXYZ = h5read(upperFile, '/geometry/vertices');
time = h5read(upperFile, '/time');
time = reshape(time, [1, size(time, 3)]);

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
NXYZs_upper = nodalXYZ(1:2, NID);
nodeID = 0;
NIDD = zeros(1, size(NID, 2));
for i = 1:1:size(NID, 2)
    if NID(i) == 1
        nodeID = nodeID + 1;
        NIDD(i) = nodeID;
    end
end

Eles_upper = Connection(1:3, EID);
for i = 1:1:size(Eles_upper, 2)
    Eles_upper(:, i) = NIDD(Eles_upper(:, i))';
end

stress = stress(:, EID, :);

nOfTimeSteps = size(stress, 3);
nOfCells = size(stress, 2);
shearStress_upper = zeros(nOfCells, nOfTimeSteps);
n = normal_dir;
t = [-n(2), n(1), 0.];

% Calculate the shear stress as <n \cdot \sigma \cdot t>
for j = 1:1:nOfTimeSteps
    for i = 1:1:nOfCells
        sigma = [stress(1, i, j), stress(2, i, j), stress(3, i, j); ...
                 stress(2, i, j), stress(4, i, j), stress(5, i, j); ...
                 stress(3, i, j), stress(5, i, j), stress(6, i, j)];
        shearStress_upper(i, j) = n * sigma * t';
    end
end

%% Plot surface cell view of shear stress
videoflag = true;
framerate = 8;
figNo = 1;
timeWindow = [0, 150];
fontsize = 25;
if videoflag == true
    fig = figure(figNo);
    fig.Position = [1000, 597, 2 * fig.Position(3), 2 * fig.Position(4)];
    xrange = [-110, 110];
    yrange = [-110, 110];
    
    % Initialize names
    videoname = strcat(videoprefix, '_shearStress.mp4');
    
    % Initialize video
    myVideo = VideoWriter(strcat('Videos/', videoname), 'MPEG-4');
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
        p_upper = patch('Faces', Eles_upper' ,'Vertices', 1e3 * NXYZs_upper' ,'FaceVertexCData',shearStress_upper(:, i)/1.0e6,'FaceColor','flat');
        hold on;
        p_lower = patch('Faces', Eles_lower' ,'Vertices', 1e3 * NXYZs_lower' ,'FaceVertexCData',shearStress_lower(:, i)/1.0e6,'FaceColor','flat');
        % p = patch('XData', NXYZs(1, :), 'YData', NXYZs(2, :), 'edgecolor', 'none');
        p_upper.EdgeColor = 'none';
        p_lower.EdgeColor = 'none';
        % p.Faces = Eles';
        % p.FaceColor = 'none';
        % p.FaceVertexCData = shearStress(:, i);
        c = colorbar;
        caxis([0 8]);
        ylabel(c,'Shear Stress Along Fault-Parallel Direction [MPa]','FontName','Avenir','FontSize',fontsize);
        axis equal;
        grid on;
        xlim(xrange);
        ylim(yrange);
        xlabel('X [mm]');
        ylabel('Y [mm]');
        title(strcat('Time = ', num2str(1e6 * time(i), '%2f'), ' [\mu s]'));
        hold on;
        
        % FaultLine, wire position and VS region
        plot(1e3 * [faultST(1), faultND(1)], ...
             1e3 * [faultST(2), faultND(2)], ...
             'w', 'linewidth', 1.0);
        hold on;
        
        % Plot the VS region
        plot(1e3 * [VSstart(1), VSend(1)], ...
             1e3 * [VSstart(2), VSend(2)], ...
             'r', 'linewidth', 1.5);
        text(30, 10, VSLabel, 'color', 'red', 'Fontsize', fontsize);
        hold on;
        
        % Scatter the wire position
        scatter(1e3 * WirePos1(1), 1e3 * WirePos1(2), 80, 'filled', 'w');
        text(1e3 * WirePos1(1) + 5, 1e3 * WirePos1(2) - 2, wireLabel, 'color', 'w', 'Fontsize', fontsize);
        
        % scatter(1e3 * targetCoords(1), 1e3 * targetCoords(2), 80, 'filled', 'y');
        % text(1e3 * targetCoords(1) + 5, 1e3 * targetCoords(2) - 2, probeLabel, 'color', 'y', 'Fontsize', fontsize);
        hold off;
        
        set(gca, 'FontSize', fontsize);
        
        % Write the video
        frame = getframe(gcf);
        writeVideo(myVideo, frame);
    end
    close(myVideo);
end
