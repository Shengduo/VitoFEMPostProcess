% Script to generate fake measurement of stress, strain
clc,clear;
close all;

%% Some pre-calculated geometry info
videoprefix = '1WithWallDRS1.5_1.5ModA0.008Load5_Vw2_fw0.1_theta0.036_8000000000000000.0_NULoad2dir0';
dataFile = strcat('../dumpFiles/', videoprefix, '-domain.h5');
lowerFile = strcat('../dumpFiles/', videoprefix, '-lower_crust.h5');

% Whether or not apply average symmetry
symmetryFlag = true;

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
EID = false(1, size(Connection, 2));
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

%% Process the lower .hdf5 file for nodes in lower half
h5disp(lowerFile);
% stress = h5read(lowerFile, '/cell_fields/stress');
lower_Connection = h5read(lowerFile, '/topology/cells');
lower_Connection = lower_Connection + 1;
lower_nodalXYZ = h5read(lowerFile, '/geometry/vertices');
% time = h5read(lowerFile, '/time');
% time = reshape(time, [1, size(time, 3)]);

% Only extract the surface nodes
surface = 0.005;
lower_NID = (abs(lower_nodalXYZ(3, :) - surface) <= 1e-7);
lower_EID = false(1, size(lower_Connection, 2));
for i = 1:1:size(lower_EID, 2)
    lower_EID(i) = (sum(lower_NID(lower_Connection(:, i)')) == 3);
    if (lower_EID(i) == 1)
        lower_Connection(1:3, i) = lower_Connection(lower_NID(lower_Connection(:, i)) == 1, i);
    end
end

% Update nodes and collection
lower_NXYZs_lower = lower_nodalXYZ(1:2, lower_NID);
lower_nodeID = 0;
lower_NIDD = zeros(1, size(lower_NID, 2));
for i = 1:1:size(lower_NID, 2)
    if lower_NID(i) == 1
        lower_nodeID = lower_nodeID + 1;
        lower_NIDD(i) = lower_nodeID;
    end
end

lower_Eles_lower = lower_Connection(1:3, lower_EID);
for i = 1:1:size(lower_Eles_lower, 2)
    lower_Eles_lower(:, i) = lower_NIDD(lower_Eles_lower(:, i))';
end

%% Assign upper and lower nodes
% Assign upper nodes to be 1, lower nodes to be 0
UpperLowerID = true(1, size(NXYZs, 2));
for i = 1:1:size(lower_Eles_lower, 2)
    UpperLowerID(lower_Eles_lower(:, i)) = 0;
end
disp(strcat('Upper nodes count: ', num2str(sum(UpperLowerID))));
disp(strcat('Lower nodes count: ', num2str(size(UpperLowerID, 2) - sum(UpperLowerID))));

%% Rotation

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

% Only input the values of the observation
windowID = ((XYZ(1, :) >= VS_start_x - 0.001) & (XYZ(1, :) <= VS_end_x + 0.001) ...
            & (XYZ(2, :) >= -VS_end_y - 0.001) & (XYZ(2, :) <= VS_end_y + 0.001));
UpperID = UpperLowerID & windowID;
LowerID = (~UpperLowerID) & windowID;

DICdisp_up = zeros(2, size(X_up, 1), size(X_up, 2), size(time, 2));
DICdisp_low = zeros(2, size(X_low, 1), size(X_low, 2), size(time, 2));

% Get the values of displacements at [x, y]s
for i = 1:1:nOfTimeSteps
    % Symmetrize the displacements
    % Print timestep
    disp(strcat("Calculating displacement, time step: ", num2str(i)));
    F_displacement_up_x = scatteredInterpolant(XYZ(:, UpperID)', squeeze(displacement(1, UpperID, i))', 'linear');
    F_displacement_up_y = scatteredInterpolant(XYZ(:, UpperID)', squeeze(displacement(2, UpperID, i))', 'linear');
    F_displacement_low_x = scatteredInterpolant(XYZ(:, LowerID)', squeeze(displacement(1, LowerID, i))', 'linear');
    F_displacement_low_y = scatteredInterpolant(XYZ(:, LowerID)', squeeze(displacement(2, LowerID, i))', 'linear');
    
    DICdisp_up(1, :, :, i) = F_displacement_up_x(X_up, Y_up);
    DICdisp_up(2, :, :, i) = F_displacement_up_y(X_up, Y_up);
    DICdisp_low(1, :, :, i) = F_displacement_low_x(X_low, Y_low);
    DICdisp_low(2, :, :, i) = F_displacement_low_y(X_low, Y_low);
    
    % Apply averaging between upper and lower layer
    if symmetryFlag == true
        DICdisp_up(1, :, :, i) = (DICdisp_up(1, :, :, i) - DICdisp_low(1, :, :, i)) ./ 2;
        DICdisp_low(1, :, :, i) = - DICdisp_up(1, :, :, i);
        DICdisp_up(2, :, :, i) = (DICdisp_up(2, :, :, i) + DICdisp_low(2, :, :, i)) ./ 2;
        DICdisp_low(2, :, :, i) = DICdisp_up(2, :, :, i);
    end
end

%% Calculate strain and stress by central difference, in plane-stress case
DICstrain_up = zeros(3, size(X_up, 1), size(X_up, 2), size(time, 2));
DICstrain_low = zeros(3, size(X_low, 1), size(X_low, 2), size(time, 2));
tempstrain_up = zeros(2, size(X_up, 1), size(X_up, 2));
tempstrain_low = zeros(2, size(X_low, 1), size(X_low, 2));
DICstress_up = zeros(3, size(X_up, 1), size(X_up, 2), size(time, 2));
DICstress_low = zeros(3, size(X_low, 1), size(X_low, 2), size(time, 2));

% Get fault shear stress and normal stress
faultStress = zeros(3, size(x_up, 2), size(time, 2));
faultSlip = zeros(size(x_up, 2), size(time, 2));

% Young's modulus and Poisson's ratio
Ed = 5300e6;
vi = 0.35;

% Calculate strain
for t = 1:1:nOfTimeSteps
    % Print timestep
    disp(strcat("Calculating strain and stress, time step: ", num2str(t)));
    
    % Calculate fault slip
    faultSlip(:, t) = DICdisp_up(1, 1, :, t) - DICdisp_low(1, 1, :, t);
    
    % eps_xx and temp eps_yx
    for ii = 2:1:size(X_up, 2) - 1
        DICstrain_up(1, :, ii, t) = (DICdisp_up(1, :, ii + 1, t) - DICdisp_up(1, :, ii - 1, t)) ./ (2.0 * stepsize * pxsize);
        DICstrain_low(1, :, ii, t) = (DICdisp_low(1, :, ii + 1, t) - DICdisp_low(1, :, ii - 1, t)) ./ (2.0 * stepsize * pxsize);
        tempstrain_up(1, :, ii) = (DICdisp_up(2, :, ii + 1, t) - DICdisp_up(2, :, ii - 1, t)) ./ (2.0 * stepsize * pxsize);
        tempstrain_low(1, :, ii) = (DICdisp_low(2, :, ii + 1, t) - DICdisp_low(2, :, ii - 1, t)) ./ (2.0 * stepsize * pxsize);
    end
    DICstrain_up(1, :, 1, t) = -(3 * DICdisp_up(1, :, 1, t) - 4 * DICdisp_up(1, :, 2, t) + DICdisp_up(1, :, 3, t)) ./ (2.0 * stepsize * pxsize);
    DICstrain_up(1, :, size(X_up, 2), t) = (3 * DICdisp_up(1, :, size(X_up, 2), t) - 4 * DICdisp_up(1, :, size(X_up, 2) - 1, t) + DICdisp_up(1, :, size(X_up, 2) - 2, t)) ./ (2.0 * stepsize * pxsize);
    DICstrain_low(1, :, 1, t) = -(3 * DICdisp_low(1, :, 1, t) - 4 * DICdisp_low(1, :, 2, t) + DICdisp_low(1, :, 3, t)) ./ (2.0 * stepsize * pxsize);
    DICstrain_low(1, :, size(X_up, 2), t) = (3 * DICdisp_low(1, :, size(X_up, 2), t) - 4 * DICdisp_low(1, :, size(X_up, 2) - 1, t) + DICdisp_low(1, :, size(X_up, 2) - 2, t)) ./ (2.0 * stepsize * pxsize);
    
    tempstrain_up(1, :, 1) = -(3 * DICdisp_up(2, :, 1, t) - 4 * DICdisp_up(2, :, 2, t) + DICdisp_up(2, :, 3, t)) ./ (2.0 * stepsize * pxsize);
    tempstrain_up(1, :, size(X_up, 2)) = (3 * DICdisp_up(2, :, size(X_up, 2), t) - 4 * DICdisp_up(2, :, size(X_up, 2) - 1, t) + DICdisp_up(2, :, size(X_up, 2) - 2, t)) ./ (2.0 * stepsize * pxsize);
    tempstrain_low(1, :, 1) = -(3 * DICdisp_low(2, :, 1, t) - 4 * DICdisp_low(2, :, 2, t) + DICdisp_low(2, :, 3, t)) ./ (2.0 * stepsize * pxsize);
    tempstrain_low(1, :, size(X_up, 2)) = (3 * DICdisp_low(2, :, size(X_up, 2), t) - 4 * DICdisp_low(2, :, size(X_up, 2) - 1, t) + DICdisp_low(2, :, size(X_up, 2) - 2, t)) ./ (2.0 * stepsize * pxsize);
    
    % eps_yy and temp eps_xy
    for jj = 2:1:size(X_up, 1) - 1
        DICstrain_up(2, jj, :, t) = (DICdisp_up(2, jj + 1, :, t) - DICdisp_up(2, jj - 1, :, t)) ./ (2.0 * stepsize * pxsize);
        DICstrain_low(2, jj, :, t) = (DICdisp_low(2, jj + 1, :, t) - DICdisp_low(2, jj - 1, :, t)) ./ (2.0 * stepsize * pxsize);
        tempstrain_up(2, jj, :) = (DICdisp_up(1, jj + 1, :, t) - DICdisp_up(1, jj - 1, :, t)) ./ (2.0 * stepsize * pxsize);
        tempstrain_low(2, jj, :) = (DICdisp_low(1, jj + 1, :, t) - DICdisp_low(1, jj - 1, :, t)) ./ (2.0 * stepsize * pxsize);
    end
    DICstrain_up(2, 1, :, t) = -(3 * DICdisp_up(2, 1, :, t) - 4 * DICdisp_up(2, 2, :, t) + DICdisp_up(2, 3, :, t)) ./ (2.0 * stepsize * pxsize);
    DICstrain_up(2, size(X_up, 1), :, t) = (3 * DICdisp_up(2, size(X_up, 1), :, t) - 4 * DICdisp_up(2, size(X_up, 1) - 1, :, t) + DICdisp_up(2, size(X_up, 1) - 2, :, t)) ./ (2.0 * stepsize * pxsize);
    DICstrain_low(2, 1, :, t) = -(3 * DICdisp_low(2, 1, :, t) - 4 * DICdisp_low(2, 2, :, t) + DICdisp_low(2, 3, :, t)) ./ (2.0 * stepsize * pxsize);
    DICstrain_low(2, size(X_up, 1), :, t) = (3 * DICdisp_low(2, size(X_up, 1), :, t) - 4 * DICdisp_low(2, size(X_up, 1) - 1, :, t) + DICdisp_low(2, size(X_up, 1) - 2, :, t)) ./ (2.0 * stepsize * pxsize);
    
    tempstrain_up(2, 1, :) = -(3 * DICdisp_up(1, 1, :, t) - 4 * DICdisp_up(1, 2, :, t) + DICdisp_up(1, 3, :, t)) ./ (2.0 * stepsize * pxsize);
    tempstrain_up(2, size(X_up, 1), :) = (3 * DICdisp_up(1, size(X_up, 1), :, t) - 4 * DICdisp_up(1, size(X_up, 1) - 1, :, t) + DICdisp_up(1, size(X_up, 1) - 2, :, t)) ./ (2.0 * stepsize * pxsize);
    tempstrain_low(2, 1, :) = -(3 * DICdisp_low(1, 1, :, t) - 4 * DICdisp_low(1, 2, :, t) + DICdisp_low(1, 3, :, t)) ./ (2.0 * stepsize * pxsize);
    tempstrain_low(2, size(X_up, 1), :) = (3 * DICdisp_low(1, size(X_up, 1), :, t) - 4 * DICdisp_low(1, size(X_up, 1) - 1, :, t) + DICdisp_low(1, size(X_up, 1) - 2, :, t)) ./ (2.0 * stepsize * pxsize);
    
    % Recover eps_xy = eps_yx
    DICstrain_up(3, :, :, t) = (tempstrain_up(1, :, :) + tempstrain_up(2, :, :)) ./ 2.0; 
    DICstrain_low(3, :, :, t) = (tempstrain_low(1, :, :) + tempstrain_low(2, :, :)) ./ 2.0; 
    
    % Calculate stress
    DICstress_up(1, :, :, t) = Ed / (1 - vi^2) * (DICstrain_up(1, :, :, t) + vi * DICstrain_up(2, :, :, t));
    DICstress_up(2, :, :, t) = Ed / (1 - vi^2) * (DICstrain_up(2, :, :, t) + vi * DICstrain_up(1, :, :, t));
    DICstress_up(3, :, :, t) = Ed / (1 + vi) * (DICstrain_up(3, :, :, t));
    DICstress_low(1, :, :, t) = Ed / (1 - vi^2) * (DICstrain_low(1, :, :, t) + vi * DICstrain_low(2, :, :, t));
    DICstress_low(2, :, :, t) = Ed / (1 - vi^2) * (DICstrain_low(2, :, :, t) + vi * DICstrain_low(1, :, :, t));
    DICstress_low(3, :, :, t) = Ed / (1 + vi) * (DICstrain_low(3, :, :, t));
    
    % Calculate fault stress;
    faultStress(:, :, t) = (DICstress_up(:, 1, :, t) + DICstress_low(:, 1, :, t)) ./ 2;
end

%% Plot DIC fake fault shear stress versus time
fontsize = 25;
figNo = 1;
load('BRColorScale.mat');

% Several wave speed to show on the aplot
cp = 2662.4;
cs = 1279;
nu = 0.35;
cr = (0.874 + 0.196 * nu - 0.043 * nu^2 - 0.055 * nu^3) * cs;
cX = [60, 80];
crY = [10, (cX(2) - cX(1)) * 1e3 / cr + 10];
csY = [10, (cX(2) - cX(1)) * 1e3 / cs + 10];
cpY = [10, (cX(2) - cX(1)) * 1e3 / cp + 10];

% Adjust to compare with plots from simulation results
plotTime = time - 10e-6;
nOfTimeSteps = size(plotTime, 2);
plot_x_up = x_up + norm(WirePos1, 2);
VSregion = 1e3 * [norm(VSstart - WirePos1, 2), norm(VSend - WirePos1, 2)];
VitoColorFlag = 1;

% Prestress
P = 14.3e6;
alpha = 29 * pi / 180;
sigma_Pre = -P * [sin(alpha)^2; cos(alpha)^2; sin(alpha) * cos(alpha)];
totalFaultStress = faultStress + sigma_Pre;
%% Save a X-T diagram plot of shear stress (only observing window)
plotflag = true;
if plotflag == true
    fig = figure(figNo);
    % Trange = [0, 150];
    Trange = [30, 110];
    Xrange = [VSregion(1), VSregion(1) + 45];
    fig.Position(3:4) = 1.5 * fig.Position(3:4);
    
    % Initialize names
    plotname = strcat(pwd, '/../plots/', videoprefix, '_X-TofShearstress_window_DIC.png');
    
    % Plot sliprate on X-T
    [Tsteps, Xsteps] = meshgrid(1e6 * plotTime, 1e3 * plot_x_up);
    %subplot(2,2,iii)
    h = pcolor(Xsteps', Tsteps', (squeeze(- totalFaultStress(3, :, :)) ./ 1e6)');
    shading interp;
    if VitoColorFlag == 1
        plotname = strcat(pwd, '/../Vitoplots/', videoprefix, '_X-TofShearStress_window_DIC.png');
        colormap(black_rainbow_plus_long);
    end
    hold on;
    xline(VSregion(1), 'r' ,'linewidth', 2.0);
    xline(VSregion(2), 'r' ,'linewidth', 2.0);
    text(VSregion(1)+ 5, 40, 'VS region', 'color', 'r', 'Fontsize', fontsize);
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
    title('X-T diagram of Shear stress');
    set(gca, 'FontSize', fontsize);
    
    % Save the figure
    print(figure(figNo) ,plotname, '-dpng', '-r500');
end
figNo = figNo + 1;

