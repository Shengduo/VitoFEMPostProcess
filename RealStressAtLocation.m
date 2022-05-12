% Compute real stress at distance above and into the gouge
clc,clear;
close all;
stress_distance = 3;
z_location = 0.0049;

%% Some pre-calculated geometry info
videoprefix = '1WithWallDRS1.5_1.5ModA0.008AmB0.005Load5_Vw2_fw0.1_theta0.036_-11_NULoad2dir0';
upperFile = strcat('../dumpFiles/', videoprefix, '-upper_crust.h5');

% Whether or not apply average symmetry
symmetryFlag = false;

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
VSregion = 1e3 * [norm(VSstart - WirePos1, 2), norm(VSend - WirePos1, 2)];

% Several wave speed to show on the aplot
cp = 2662.4;
cs = 1279;
nu = 0.35;
cr = (0.874 + 0.196 * nu - 0.043 * nu^2 - 0.055 * nu^3) * cs;
cX = [60, 80];
crY = [10, (cX(2) - cX(1)) * 1e3 / cr + 10];
csY = [10, (cX(2) - cX(1)) * 1e3 / cs + 10];
cpY = [10, (cX(2) - cX(1)) * 1e3 / cp + 10];

%% Process the domain file for displacements
h5disp(upperFile);
time = h5read(upperFile, '/time');
time = reshape(time, [1, size(time, 3)]);
stress = h5read(upperFile, '/cell_fields/stress');
connection = h5read(upperFile, '/topology/cells');
connection = connection + 1;
nodalXYZ = h5read(upperFile, '/geometry/vertices');


%% Rotation
% Apply rotation to displacement and nodal XYZ to align with the 
% fault coordinate system
Q = [tangential_dir(1:3); normal_dir(1:3); [0, 0, 1]];

% Rotate nodal xyzs
XYZ = Q * nodalXYZ;

% Get the mesh stored as triangulation
TR = triangulation(connection', XYZ');

%% Get the interpolation of slip at these locations
% Get the start and end location of  VS region
VS_start_x = norm(VSstart - 0, 2);
VS_end_x = norm(VSend - 0, 2);

% Camera specific values
Pix_10mm = 81;  %Number of pixels in a grid of 10 mm ;  144 pixels in 2 mm -> 72 pixels in 1 mm
pxsize=10/(Pix_10mm)*1e-3; % pixel size in m/pixel  (Pix_10mm pixels  in 10 mm)
stepsize = 1;

VS_start_y = pxsize;
VS_end_y = 20e-3;

% First specific the x and y's
x_up = VS_start_x : stepsize * pxsize : VS_end_x;
x_low = x_up;
y_up = VS_start_y : stepsize * pxsize : VS_end_y;     % Assume the window height is 20 mm
y_low = - y_up;

[X_up, Y_up] = meshgrid(x_up, y_up);
[X_low, Y_low] = meshgrid(x_low, y_low);

% Get the querying positions
[Xq, Yq] = meshgrid(x_up, stress_distance / 1e3);
P = [Xq', Yq', z_location * ones(size(Xq, 2), 1)];

% Get which elements the points are at
elementID = pointLocation(TR, P);

% Allocate memory for results
totalStressAtDist = zeros(size(elementID, 1), size(time, 2));

% Initialize figNo
figNo = 1;
load('BRColorScale.mat');
fontsize = 25;

% Adjust to compare with plots from simulation results
plotTime = time - 10e-6;
nOfTimeSteps = size(plotTime, 2);
plot_x_up = x_up + norm(WirePos1, 2);
VSregion = 1e3 * [norm(VSstart - WirePos1, 2), norm(VSend - WirePos1, 2)];
VitoColorFlag = 1;

%% Save a X-T diagram plot of shear stress above the surface at stress_distance (mm)
plotflag = false;
if plotflag == true
    fig = figure(figNo);
    % Trange = [0, 150];
    Trange = [30, 110];
    Xrange = [VSregion(1), VSregion(1) + 45];
    fig.Position(3:4) = 1.5 * fig.Position(3:4);

    % Initialize names
    plotname = strcat(pwd, '/../plots/', videoprefix, '_X-TofShearstress_window_', num2str(stress_distance), '_', num2str(z_location), '.png');
    stress_at_dist = zeros(size(x_up, 2), size(time, 2));
    for t = 1:1:size(time, 2)
        for ele = 1:1:size(elementID, 1)
            tempvecstress = stress(:, elementID(ele), t);
            tempmatstress = tempvecstress([1, 6, 5; 6, 2, 4; 5, 4, 3]);
            tempmatstress = Q * tempmatstress * Q';
            totalStressAtDist(ele, t) = tempmatstress(1, 2);
        end        
    end
    % Plot sliprate on X-T
    [Tsteps, Xsteps] = meshgrid(1e6 * plotTime, 1e3 * plot_x_up);
    %subplot(2,2,iii)
    h = pcolor(Xsteps', Tsteps', (squeeze(- totalStressAtDist) ./ 1e6)');
    shading interp;
    if VitoColorFlag == 1
        plotname = strcat(pwd, '/../Vitoplots/', videoprefix, '_X-TofShearstress_window_', num2str(stress_distance), '_', num2str(z_location), '.png');
        colormap(black_rainbow_plus_long);
    end
    hold on;
    xline(VSregion(1), 'r' ,'linewidth', 2.0);
    xline(VSregion(2), 'r' ,'linewidth', 2.0);
    text(VSregion(1)+ 5, 40, 'Gouge Region', 'color', 'r', 'Fontsize', fontsize);
    % Add the wave speeds

    cX = [55, 65];
    crY = [40, (cX(2) - cX(1)) * 1e3 / cr + 40];
    csY = [40, (cX(2) - cX(1)) * 1e3 / cs + 40];
    cpY = [40, (cX(2) - cX(1)) * 1e3 / cp + 40];

    plot(cX, crY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, crY(2)+2, strcat('$c_r$ = 1.20 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 10, 'interpreter', 'latex');
    plot(cX, csY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, csY(2) - 1, strcat('$c_s$ = 1.28 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 10, 'interpreter', 'latex');
    plot(cX, cpY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, cpY(2), strcat('$c_p$ = 2.66 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 10, 'interpreter', 'latex');
    hold off;
    set(h, 'EdgeColor', 'None');
    c = colorbar;
    caxis([2, 10]);
    ylabel(c,'Shear stress [MPa]','FontName','Avenir','FontSize',fontsize);
    xlim(Xrange);
    ylim(Trange);
    xlabel('Distance along the fault [mm]');
    ylabel('Time [\mus]');
    title(['X-T Shear stress at ', num2str(stress_distance), ' ', 'mm']);
    set(gca, 'FontSize', fontsize);

    % Save the figure
    print(figure(figNo) ,plotname, '-dpng', '-r500');
end
figNo = figNo + 1;

%% Front surf file to plot displacement
% Get the displacement at these locations
domainFile = strcat('../dumpFiles/', videoprefix, '-domain.h5');
h5disp(domainFile);

% Get lower element ID's
lowerFile = strcat('../dumpFiles/', videoprefix, '-lower_crust.h5');
connection_low = h5read(lowerFile, '/topology/cells');
connection_low = connection_low + 1;
TR_low = triangulation(connection_low', XYZ');
P_low = [Xq', -Yq', z_location * ones(size(Xq, 2), 1)];
elementID_low = pointLocation(TR_low, P_low);

% Pre-allocate memory for the results
faultSlipRate = zeros(size(elementID, 1), size(time, 2));
velocity = h5read(domainFile, '/vertex_fields/velocity');
for t = 1:1:size(time, 2)
    velocity(:, :, t) = Q * velocity(:, :, t);
end
%% Save a X-T diagram plot of surface slip rate at stress_distance, z_location (observing window, smaller range)
plotflag = true;
if plotflag == true
    fig = figure(figNo);
    % Trange = [0, 150];
    Trange = [30, 110];
    Xrange = [VSregion(1), VSregion(1) + 45];
    fig.Position(3:4) = 1.5 * fig.Position(3:4);
    
    % Compute the slip rate
    for t = 1:1:size(time, 2)
        for ele = 1:1:size(elementID, 1)
            F_vel_up_x = scatteredInterpolant(XYZ(:, connection(:, elementID(ele)))', squeeze(velocity(1, connection(:, elementID(ele)), t))', 'natural');
            F_vel_low_x = scatteredInterpolant(XYZ(:, connection_low(:, elementID_low(ele)))', squeeze(velocity(1, connection_low(:, elementID_low(ele)), t))', 'natural');
            faultSlipRate(ele, t) = F_vel_up_x(x_up(ele), stress_distance / 1e3, z_location) - F_vel_low_x(x_up(ele), -stress_distance / 1e3, z_location);
        end
    end
    
    % Initialize names
    plotname = strcat(pwd, '/../plots/', videoprefix, '_X-TofSlipRate_window_small_', num2str(stress_distance), '_', num2str(z_location), '.png');
    
    % Plot sliprate on X-T
    [Tsteps, Xsteps] = meshgrid(1e6 * plotTime, 1e3 * plot_x_up);
    h = pcolor(Xsteps', Tsteps', (-faultSlipRate)');
    shading interp;
    
    if VitoColorFlag == 1
        plotname = strcat(pwd, '/../Vitoplots/', videoprefix, '_X-TofSlipRate_window_small_', num2str(stress_distance), '_', num2str(z_location), '.png');
        colormap(flipud(black_rainbow_shear_long));
    end
    hold on;
    xline(VSregion(1), 'r' ,'linewidth', 2.0);
    xline(VSregion(2), 'r' ,'linewidth', 2.0);
    text(VSregion(1)+ 5, 40, 'VS region', 'color', 'r', 'Fontsize', fontsize);
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
    xlabel('Distance along the fault [mm]');
    ylabel('Time [\mus]');
    title(['X-T diagram of Slip rate at ', ' ', num2str(stress_distance), ' ', '[mm]']);
    set(gca, 'FontSize', fontsize);
    
    % Save the figure
    print(figure(figNo) ,plotname, '-dpng', '-r500');
end
figNo = figNo + 1;
