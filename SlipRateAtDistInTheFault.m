clc,clear;
close all;
% This script process the faultfile and get the slip at certain location

load('BRColorScale.mat');
VitoColorFlag = true;
tractionOffsetFlag = false;

%% Input parameters
videoprefix = '1WithWallDRS1.5_1.5ModA0.008AmB0.005Load5_Vw2_fw0.1_theta0.036_-11_NULoad2dir0';
% Select the surface nodes
Distance_To_Surface = 0.005;

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

%% Get nodal information
nOfNodes = size(nodalXYZ, 2);

% Modify videoprefix, show number of nodes
videoprefix = strcat(videoprefix, "_Distance", num2str(Distance_To_Surface));

% Read nodal slip, slip rate and traction
SlipRate = h5read(faultFileName, '/vertex_fields/slip_rate');
connection = h5read(faultFileName, '/topology/cells');
connection = connection + 1;

% Read nodal slip
% Input the first wire position
WirePos1 = [-0.025657; -0.014222; 0];
VSstart = [0.006354, 0.003522, 0]';
VSend = [0.063204, 0.035034, 0]';
VSregion = 1e3 * [norm(VSstart - WirePos1, 2), norm(VSend - WirePos1, 2)];
% VSregion = [50, 120];

% Calculate NodalXYZ2D 
NodalXYZ2D = zeros(2, size(nodalXYZ, 2));
for i = 1:1:size(nodalXYZ, 2)
    NodalXYZ2D(1, i) = sign(nodalXYZ(1, i) - WirePos1(1)) * norm(nodalXYZ(1:2, i) - WirePos1(1:2), 2);
    NodalXYZ2D(2, i) = nodalXYZ(3, i);
end

% Fault Range
xrange = [-100, 150];
Vrange = [0, 20];

% Magnitude of slip rate
slipRateMag = zeros(nOfNodes, nOfTimeSteps);

for t = 1:1:nOfTimeSteps
    for i = 1:1:nOfNodes
        slipRateMag(i, t) = norm(SlipRate(:, i, t), 2);
    end
end
figNo = 1;

% Get the mesh stored as triangulation
TR = triangulation(connection', NodalXYZ2D');

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

% Get the querying positions
[Xq, Yq] = meshgrid(x_up, Distance_To_Surface);
P = [Xq', Yq'];

elementID = pointLocation(TR, P);
SlipRateAtDist = zeros(size(Xq, 2), size(time, 2));

for t = 1:1:size(time, 2)
    for ele = 1:1:size(elementID, 1)
        vel_x = scatteredInterpolant(NodalXYZ2D(:, connection(:, elementID(ele)))', squeeze(slipRateMag(connection(:, elementID(ele)), t)), 'natural');
        SlipRateAtDist(ele, t) = vel_x(x_up(ele), Distance_To_Surface);
    end
end


%% Save a X-T diagram plot of surface slip rate at stress_distance, z_location (observing window, smaller range)
plotflag = true;
if plotflag == true
    fig = figure(figNo);
    % Trange = [0, 150];
    Trange = [30, 110];
    Xrange = [VSregion(1), VSregion(1) + 45];
    fig.Position(3:4) = 1.5 * fig.Position(3:4);

    % Initialize names
    plotname = strcat(pwd, '/../plots/', videoprefix, '_X-TofSlipRate_window_small_surface_', num2str(Distance_To_Surface), '.png');

    % Plot sliprate on X-T
    [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * x_up);
    h = pcolor(Xsteps', Tsteps', SlipRateAtDist');
    shading interp;

    if VitoColorFlag == 1
        plotname = strcat(pwd, '/../Vitoplots/', videoprefix, '_X-TofSlipRate_window_small_surface_', num2str(Distance_To_Surface), '.png');
        colormap(flipud(black_rainbow_shear_long));
    end
    hold on;
    xline(VSregion(1), 'r' ,'linewidth', 2.0);
    xline(VSregion(2), 'r' ,'linewidth', 2.0);
    text(VSregion(1)+ 5, 40, 'VS region', 'color', 'r', 'Fontsize', fontsize, 'interpreter', 'latex');
    % Add the wave speeds

    cX = [55, 65];
    crY = [60, (cX(2) - cX(1)) * 1e3 / cr + 60];
    csY = [60, (cX(2) - cX(1)) * 1e3 / cs + 60];
    cpY = [60, (cX(2) - cX(1)) * 1e3 / cp + 60];

    plot(cX, crY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, crY(2)+2, strcat('$c_r$ = 1.20 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 10, 'interpreter', 'latex');
    plot(cX, csY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, csY(2) - 1, strcat('$c_s$ = 1.28 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 10, 'interpreter', 'latex');
    plot(cX, cpY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, cpY(2), strcat('$c_p$ = 2.66 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 10, 'interpreter', 'latex');
    hold off;
    set(h, 'EdgeColor', 'None');
    c = colorbar;
    caxis([0, 0.8]);
    ylabel(c,'Slip rate [m/s]','FontName','Avenir','FontSize',fontsize);
    xlim(Xrange);
    ylim(Trange);
    xlabel('Distance along the fault [mm]', 'interpreter', 'latex');
    ylabel('Time [$\mu$s]', 'interpreter', 'latex');
    title(['X-T diagram of Slip rate at Z = ', ' ', num2str(1e3 * Distance_To_Surface, '%.0f'), ' ', '[mm]'], 'interpreter', 'latex');
    set(gca, 'FontSize', fontsize);

    % Save the figure
    print(figure(figNo) ,plotname, '-dpng', '-r500');
end
figNo = figNo + 1;


