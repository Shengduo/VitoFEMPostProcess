% Set environments
setEnvironment;
videoprefix = 'W8_4ExcCorr2NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2e+16_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8';
% videoprefix = 'W8_4ExcCorr2NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load8_Vw2_fw0.1_theta0.8061_-9_NULoad2dir0_duration200_8';
% videoprefix = "W12Hom0.0015_0.47_SW_0.1_2_8_36_40_W8_4ExcCorr2NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load8_Vw2_fw0.1_theta0.2208_-9_NULoad2dir0_duration200_50";
% videoprefix = "BlkW8_SW_0.01_3.56_4_8_36_40_7e-7_W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.012AmB0.01Load8_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_12";
upperFile = strcat('../dumpFiles/', videoprefix, '-upper_crust.h5');
lowerFile = strcat('../dumpFiles/', videoprefix, '-lower_crust.h5');
frontFileName = strcat('../frontsurfFiles/', videoprefix, '-frontsurf.h5');

load('BRColorScale.mat');

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
time = time - 10e-6;

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
normalStress_lower = zeros(nOfCells, nOfTimeSteps);
n = normal_dir;
t = [-n(2), n(1), 0.];

% Calculate the shear stress as <n \cdot \sigma \cdot t>
for j = 1:1:nOfTimeSteps
    for i = 1:1:nOfCells
        sigma = [stress(1, i, j), stress(4, i, j), stress(6, i, j); ...
                 stress(4, i, j), stress(2, i, j), stress(5, i, j); ...
                 stress(6, i, j), stress(5, i, j), stress(3, i, j)];
        shearStress_lower(i, j) = n * sigma * t';
        normalStress_lower(i, j) = -n * sigma * n';
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
time = time - 10e-6;
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
normalStress_upper = zeros(nOfCells, nOfTimeSteps);

n = normal_dir;
t = [-n(2), n(1), 0.];

% Calculate the shear stress as <n \cdot \sigma \cdot t>
for j = 1:1:nOfTimeSteps
    for i = 1:1:nOfCells
        sigma = [stress(1, i, j), stress(4, i, j), stress(6, i, j); ...
                 stress(4, i, j), stress(2, i, j), stress(5, i, j); ...
                 stress(6, i, j), stress(5, i, j), stress(3, i, j)];
        shearStress_upper(i, j) = n * sigma * t';
        normalStress_upper(i, j) = -n * sigma * n';
    end
end

%% Plot surface cell view of shear and normal stress
videoflag = true;
framerate = 8;
figNo = 1;
timeWindow = [0, 150];
fontsize = 25;
if videoflag == true
    fig = figure(figNo);
    % set(gcf, 'Position', [10, 10, 1900, 1409])
    figPos = [1900, 1409];
    fig.Position(3:4) = figPos;
    xrange = [-100, 100];
    yrange = [-100, 100];
    hold on;
    fig.Position(3:4) = figPos;
    % Initialize names
%     videoname = strcat(videoprefix, '_shearStress.avi');
    
%     % Initialize video
%     myVideo = VideoWriter(strcat('../stress_20230221/', videoname), 'Motion JPEG AVI');
%     myVideo.FrameRate = framerate;
%     myVideo.Quality = 100;
%     open(myVideo);
    plotPrefix_shear = "../stress_20230221/" + videoprefix + "_shearStress";
    plotPrefix_normal = "../stress_20230221/" + videoprefix + "_normalStress";
    % rmdir(plotPrefix);
    
    mkdir(plotPrefix_shear);
    mkdir(plotPrefix_normal);

    %% Shoot the shear stress video
    for i = 1:1:size(time, 2)
        hold on;
        % fig.Position(3:4) = figPos;
        
        if time(i) * 1e6 > timeWindow(2)
            break;
        end
        p_upper = patch('Faces', Eles_upper' ,'Vertices', 1e3 * NXYZs_upper' ,'FaceVertexCData',shearStress_upper(:, i)/1.0e6,'FaceColor','flat');
        hold on;
        p_lower = patch('Faces', Eles_lower' ,'Vertices', 1e3 * NXYZs_lower' ,'FaceVertexCData',shearStress_lower(:, i)/1.0e6,'FaceColor','flat');
        % p = patch('XData', NXYZs(1, :), 'YData', NXYZs(2, :), 'edgecolor', 'none');
        p_upper.EdgeColor = 'none';
        p_lower.EdgeColor = 'none';
        colormap(black_rainbow_plus_long);
        % p.Faces = Eles';
        % p.FaceColor = 'none';
        % p.FaceVertexCData = shearStress(:, i);
        c = colorbar;
        mid = 14.3 * cosd(29) * sind(29);
        crange2 = 10;
        clim([crange2 - 2 * (crange2 - mid), crange2]);
        c.Ticks = [crange2 - 2 * (crange2 - 6), 6, crange2];
        ylabel(c,'Shear Stress Along Fault-Parallel Direction [MPa]','FontName','Avenir','FontSize',fontsize, 'Interpreter', 'latex');
        axis equal;
        grid on;
        xlim(xrange);
        ylim(yrange);
        xlabel('X [mm]');
        ylabel('Y [mm]');
        title("Time = " + num2str(1e6 * time(i), '%.2f') + " [$\mathrm{\mu s}$]");
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

        
        set(gca, 'FontSize', fontsize);
        set(gcf, 'color', 'w');
%         plotname = plotPrefix + "/" + num2str(i) + ".png";
%         print(figure(figNo), plotname, '-dpng', '-r100');
        % Write the video
        % fig.Position(3:4) = figPos;
        frame = getframe(gcf);
        if time(i) * 1e6 >= timeWindow(1)
            % writeVideo(myVideo, frame);
            plotname = plotPrefix_shear + "/" + num2str(i) + ".png";
            print(figure(figNo), plotname, '-dpng', '-r100');
        end

        clf(fig);
    end
    close(fig);
%     close(myVideo);
    

    %% Shoot the normal stress video
    fig = figure(figNo);
    fig.Position(3:4) = figPos;
    xrange = [-100, 100];
    yrange = [-100, 100];
    hold on;
    fig.Position(3:4) = figPos;

    for i = 1:1:size(time, 2)
        hold on;
        % fig.Position(3:4) = figPos;
        
        if time(i) * 1e6 > timeWindow(2)
            break;
        end
        p_upper = patch('Faces', Eles_upper' ,'Vertices', 1e3 * NXYZs_upper' ,'FaceVertexCData',normalStress_upper(:, i)/1.0e6,'FaceColor','flat');
        hold on;
        p_lower = patch('Faces', Eles_lower' ,'Vertices', 1e3 * NXYZs_lower' ,'FaceVertexCData',normalStress_lower(:, i)/1.0e6,'FaceColor','flat');
        % p = patch('XData', NXYZs(1, :), 'YData', NXYZs(2, :), 'edgecolor', 'none');
        p_upper.EdgeColor = 'none';
        p_lower.EdgeColor = 'none';
        colormap(black_rainbow_plus_long);
        % p.Faces = Eles';
        % p.FaceColor = 'none';
        % p.FaceVertexCData = shearStress(:, i);
        c = colorbar;
        mid = 14.3 * cosd(29) * cosd(29);
        crange2 = 20;
        clim([crange2 - 2 * (crange2 - mid), crange2]);
        c.Ticks = [crange2 - 2 * (crange2 - 11), 11, crange2];
        ylabel(c,'Shear Stress Along Fault-Normal Direction [MPa]','FontName','Avenir','FontSize',fontsize, 'Interpreter', 'latex');
        axis equal;
        grid on;
        xlim(xrange);
        ylim(yrange);
        xlabel('X [mm]');
        ylabel('Y [mm]');
        title("Time = " + num2str(1e6 * time(i), '%.2f') + " [$\mathrm{\mu s}$]");
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

        
        set(gca, 'FontSize', fontsize);
        set(gcf, 'color', 'w');
%         plotname = plotPrefix + "/" + num2str(i) + ".png";
%         print(figure(figNo), plotname, '-dpng', '-r100');
        % Write the video
        % fig.Position(3:4) = figPos;
        frame = getframe(gcf);
        if time(i) * 1e6 >= timeWindow(1)
            % writeVideo(myVideo, frame);
            plotname = plotPrefix_normal + "/" + num2str(i) + ".png";
            print(figure(figNo), plotname, '-dpng', '-r100');
        end

        clf(fig);
    end
    close(fig);
%     close(myVideo);
end


%% Plot surface cell view of shear and normal stress
videoflag = false;
framerate = 8;
figNo = 1;
timeWindow = [0, 150];
fontsize = 25;
if videoflag == true
    fig = figure(figNo);
    % set(gcf, 'Position', [10, 10, 1900, 1409])
    figPos = [1900, 1409];
    fig.Position(3:4) = figPos;
    xrange = [-100, 100];
    yrange = [-100, 100];
    hold on;
    fig.Position(3:4) = figPos;
    % Initialize names
%     videoname = strcat(videoprefix, '_shearStress.avi');
    
%     % Initialize video
%     myVideo = VideoWriter(strcat('../stress_20230221/', videoname), 'Motion JPEG AVI');
%     myVideo.FrameRate = framerate;
%     myVideo.Quality = 100;
%     open(myVideo);
    plotPrefix_shear = "../stress_20230221/" + videoprefix + "_ShearByNormalStress";
    plotPrefix_normal = "../stress_20230221/" + videoprefix + "_ShearByNormalStressChange";
    % rmdir(plotPrefix);
    
    mkdir(plotPrefix_shear);
    mkdir(plotPrefix_normal);

    %% Shoot the shear stress video
    for i = 1:1:size(time, 2)
        hold on;
        % fig.Position(3:4) = figPos;
        
        if time(i) * 1e6 > timeWindow(2)
            break;
        end
        p_upper = patch('Faces', Eles_upper' ,'Vertices', 1e3 * NXYZs_upper' ,'FaceVertexCData',shearStress_upper(:, i) ./ normalStress_upper(:, i),'FaceColor','flat');
        hold on;
        p_lower = patch('Faces', Eles_lower' ,'Vertices', 1e3 * NXYZs_lower' ,'FaceVertexCData',shearStress_lower(:, i) ./ normalStress_lower(:, i),'FaceColor','flat');
        % p = patch('XData', NXYZs(1, :), 'YData', NXYZs(2, :), 'edgecolor', 'none');
        p_upper.EdgeColor = 'none';
        p_lower.EdgeColor = 'none';
        colormap(black_rainbow_plus_long);
        % p.Faces = Eles';
        % p.FaceColor = 'none';
        % p.FaceVertexCData = shearStress(:, i);
        c = colorbar;
        mid = sind(29) / cosd(29);
        crange2 = 1.;
        clim([crange2 - 2 * (crange2 - mid), crange2]);
        c.Ticks = [crange2 - 2 * (crange2 - mid), mid, crange2];
        ylabel(c,'$\tau/\sigma$','FontName','Avenir','FontSize',fontsize, 'Interpreter', 'latex');
        axis equal;
        grid on;
        xlim(xrange);
        ylim(yrange);
        xlabel('X [mm]');
        ylabel('Y [mm]');
        title("Time = " + num2str(1e6 * time(i), '%.2f') + " [$\mathrm{\mu s}$]");
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

        
        set(gca, 'FontSize', fontsize);
        set(gcf, 'color', 'w');
%         plotname = plotPrefix + "/" + num2str(i) + ".png";
%         print(figure(figNo), plotname, '-dpng', '-r100');
        % Write the video
        % fig.Position(3:4) = figPos;
        frame = getframe(gcf);
        if time(i) * 1e6 >= timeWindow(1)
            % writeVideo(myVideo, frame);
            plotname = plotPrefix_shear + "/" + num2str(i) + ".png";
            print(figure(figNo), plotname, '-dpng', '-r100');
        end

        clf(fig);
    end
    close(fig);
%     close(myVideo);
    

%     %% Shoot the delta shear / normal stress video
%     fig = figure(figNo);
%     fig.Position(3:4) = figPos;
%     xrange = [-100, 100];
%     yrange = [-100, 100];
%     hold on;
%     fig.Position(3:4) = figPos;
% 
%     for i = 1:1:size(time, 2)
%         hold on;
%         % fig.Position(3:4) = figPos;
%         
%         if time(i) * 1e6 > timeWindow(2)
%             break;
%         end
%         p_upper = patch('Faces', Eles_upper' ,'Vertices', 1e3 * NXYZs_upper' ,'FaceVertexCData',normalStress_upper(:, i)/1.0e6,'FaceColor','flat');
%         hold on;
%         p_lower = patch('Faces', Eles_lower' ,'Vertices', 1e3 * NXYZs_lower' ,'FaceVertexCData',normalStress_lower(:, i)/1.0e6,'FaceColor','flat');
%         % p = patch('XData', NXYZs(1, :), 'YData', NXYZs(2, :), 'edgecolor', 'none');
%         p_upper.EdgeColor = 'none';
%         p_lower.EdgeColor = 'none';
%         colormap(black_rainbow_plus_long);
%         % p.Faces = Eles';
%         % p.FaceColor = 'none';
%         % p.FaceVertexCData = shearStress(:, i);
%         c = colorbar;
%         mid = 14.3 * cosd(29) * cosd(29);
%         crange2 = 20;
%         clim([crange2 - 2 * (crange2 - mid), crange2]);
%         c.Ticks = [crange2 - 2 * (crange2 - 11), 11, crange2];
%         ylabel(c,'Stress Along Fault-Normal Direction [MPa]','FontName','Avenir','FontSize',fontsize, 'Interpreter', 'latex');
%         axis equal;
%         grid on;
%         xlim(xrange);
%         ylim(yrange);
%         xlabel('X [mm]');
%         ylabel('Y [mm]');
%         title("Time = " + num2str(1e6 * time(i), '%.2f') + " [$\mathrm{\mu s}$]");
%         hold on;
%         
%         % FaultLine, wire position and VS region
%         plot(1e3 * [faultST(1), faultND(1)], ...
%              1e3 * [faultST(2), faultND(2)], ...
%              'w', 'linewidth', 1.0);
%         hold on;
%         
%         % Plot the VS region
%         plot(1e3 * [VSstart(1), VSend(1)], ...
%              1e3 * [VSstart(2), VSend(2)], ...
%              'r', 'linewidth', 1.5);
%         text(30, 10, VSLabel, 'color', 'red', 'Fontsize', fontsize);
%         hold on;
%         
%         % Scatter the wire position
%         scatter(1e3 * WirePos1(1), 1e3 * WirePos1(2), 80, 'filled', 'w');
%         text(1e3 * WirePos1(1) + 5, 1e3 * WirePos1(2) - 2, wireLabel, 'color', 'w', 'Fontsize', fontsize);
%         
%         % scatter(1e3 * targetCoords(1), 1e3 * targetCoords(2), 80, 'filled', 'y');
%         % text(1e3 * targetCoords(1) + 5, 1e3 * targetCoords(2) - 2, probeLabel, 'color', 'y', 'Fontsize', fontsize);
% 
%         
%         set(gca, 'FontSize', fontsize);
%         set(gcf, 'color', 'w');
% %         plotname = plotPrefix + "/" + num2str(i) + ".png";
% %         print(figure(figNo), plotname, '-dpng', '-r100');
%         % Write the video
%         % fig.Position(3:4) = figPos;
%         frame = getframe(gcf);
%         if time(i) * 1e6 >= timeWindow(1)
%             % writeVideo(myVideo, frame);
%             plotname = plotPrefix_normal + "/" + num2str(i) + ".png";
%             print(figure(figNo), plotname, '-dpng', '-r100');
%         end
% 
%         clf(fig);
%     end
%     close(fig);
%     close(myVideo);
end



%% Grab the velocities from frontsurface
% Read time
time = h5read(frontFileName, '/time');
time = reshape(time, [1, size(time, 3)]);
time = time - 10e-6;

% Read node geometry
nodalXYZ = h5read(frontFileName, '/geometry/vertices');
nodalXYZ = nodalXYZ(1:2, :);

% Read nodal slip, slip rate and traction
velocity = h5read(frontFileName, '/vertex_fields/velocity');
velocity = velocity(1:2, :, :);
Connection = h5read(frontFileName, '/topology/cells');

% Rotate into fault coordinate system
alpha = 29 / 180 * pi;
Q = [cos(alpha), sin(alpha); -sin(alpha), cos(alpha)];

% Rotate the velocity field
for i = 1:1:size(velocity, 3)
    velocity(:, :, i) = Q * velocity(:,:,i);
end

% Rotate nodal XYZ
QnodalXYZ = Q * nodalXYZ;

plotPrefix_shear = "../stress_20230221/" + videoprefix + "_parallelVelocity";
% plotPrefix_normal = "../stress_20230221/" + videoprefix + "_normalVelocity";
% rmdir(plotPrefix);

mkdir(plotPrefix_shear);
% mkdir(plotPrefix_normal);

%% Plot fault-parallel velocities
videoflag = false;
framerate = 8;
figNo = 1;
timeWindow = [0, 150];
fontsize = 25;

if videoflag == false
    fig = figure(figNo);
    % set(gcf, 'Position', [10, 10, 1900, 1409])
    figPos = [1900, 1409];
    fig.Position(3:4) = figPos;
    xrange = [-100, 100];
    yrange = [-100, 100];
    hold on;
    fig.Position(3:4) = figPos;
    
    for i = 1 : 1 : length(time)
        hold on;
        % fig.Position(3:4) = figPos;
        
        if time(i) * 1e6 > timeWindow(2)
            break;
        end

        p = patch(1e3 * nodalXYZ(1, :), 1e3 * nodalXYZ(2, :), ...
            velocity(1, :, i), 'EdgeColor', 'none', 'FaceColor', 'interp');
        p.Faces = Connection' + 1;
        colormap((black_rainbow_plus_long));
        
        c = colorbar;
        mid = 0;
        crange2 = 7.5;
        clim([crange2 - 2 * (crange2 - mid), crange2]);
        c.Ticks = [crange2 - 2 * (crange2 - 0), 0, crange2];
        ylabel(c,'Fault-parallel velocity [m/s]','FontName','Avenir','FontSize',fontsize, 'Interpreter', 'latex');
        axis equal;
        grid on;
        xlim(xrange);
        ylim(yrange);
        xlabel('X [mm]');
        ylabel('Y [mm]');
        title("Time = " + num2str(1e6 * time(i), '%.2f') + " [$\mathrm{\mu s}$]");
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

        
        set(gca, 'FontSize', fontsize);
        set(gcf, 'color', 'w');
%         plotname = plotPrefix + "/" + num2str(i) + ".png";
%         print(figure(figNo), plotname, '-dpng', '-r100');
        % Write the video
        % fig.Position(3:4) = figPos;
        frame = getframe(gcf);
        if time(i) * 1e6 >= timeWindow(1)
            % writeVideo(myVideo, frame);
            plotname = plotPrefix_shear + "/" + num2str(i) + ".png";
            print(figure(figNo), plotname, '-dpng', '-r100');
        end

        clf(fig);
        
    end
    close(fig);
end

%% Plot fault-parallel velocities in log scale
plotPrefix_shear = "../stress_20230221/" + videoprefix + "_parallelVelocityLog";
% plotPrefix_normal = "../stress_20230221/" + videoprefix + "_normalVelocity";
% rmdir(plotPrefix);

mkdir(plotPrefix_shear);


videoflag = false;
framerate = 8;
figNo = 1;
timeWindow = [0, 150];
fontsize = 25;

if videoflag == true
    fig = figure(figNo);
    % set(gcf, 'Position', [10, 10, 1900, 1409])
    figPos = [1900, 1409];
    fig.Position(3:4) = figPos;
    xrange = [-100, 100];
    yrange = [-100, 100];
    hold on;
    fig.Position(3:4) = figPos;
    
    for i = 1 : 1 : length(time)
        hold on;
        % fig.Position(3:4) = figPos;
        
        if time(i) * 1e6 > timeWindow(2)
            break;
        end

        p = patch(1e3 * nodalXYZ(1, :), 1e3 * nodalXYZ(2, :), ...
            log10(abs(velocity(1, :, i))), 'EdgeColor', 'none', 'FaceColor', 'interp');
        p.Faces = Connection' + 1;
        % colormap((black_rainbow_plus_long));
        
        c = colorbar;
        mid = 0;
        crange2 ;
        clim([crange2 - 2 * (crange2 - mid), crange2]);
        c.Ticks = [crange2 - 2 * (crange2 - mid), mid, crange2];
        ylabel(c,'$\log(|V|)$ [m/s]','FontName','Avenir','FontSize',fontsize, 'Interpreter', 'latex');
        axis equal;
        grid on;
        xlim(xrange);
        ylim(yrange);
        xlabel('X [mm]');
        ylabel('Y [mm]');
        title("Time = " + num2str(1e6 * time(i), '%.2f') + " [$\mathrm{\mu s}$]");
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

        
        set(gca, 'FontSize', fontsize);
        set(gcf, 'color', 'w');
%         plotname = plotPrefix + "/" + num2str(i) + ".png";
%         print(figure(figNo), plotname, '-dpng', '-r100');
        % Write the video
        % fig.Position(3:4) = figPos;
        frame = getframe(gcf);
        if time(i) * 1e6 >= timeWindow(1)
            % writeVideo(myVideo, frame);
            plotname = plotPrefix_shear + "/" + num2str(i) + ".png";
            print(figure(figNo), plotname, '-dpng', '-r100');
        end

        clf(fig);
        
    end
    close(fig);
end