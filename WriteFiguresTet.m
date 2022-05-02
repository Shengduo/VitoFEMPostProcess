% Read results from hdf5 files.
clc,clear;
close all;
load('BRColorScale.mat');
VitoColorFlag = true;
tractionOffsetFlag = false;

% frontsurfFile = 'FineFHLoad15DRS1-frontsurf.h5';
% videoprefix = 'DRS1.5_8ModA0.016Load5_Vw2_fw0.1_theta0.036_NULoad2dir1';
videoprefix = '1WithWallDRS1.5_1.5ModA0.003B0.008Load5_Vw0.2_fw0.33_theta0.036_0.036_NULoad2dir0';
% videoprefix = 'ViscoElastic_theta0.043';
% videoprefix = 'DiffNULoadWithWallDRS1.5_8ModA0.016Load5_Vw2_fw0.1_theta0.036_8_NULoad2dir-1';
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

%%
% Select the surface nodes
Distance_To_Surface = 0.005;
% Distance_To_Surface = 0.004135;
% Distance_To_Surface = 0.003270;
% Distance_To_Surface = 0.002406;
I1 = (abs(nodalXYZ(3, :) - Distance_To_Surface) <= 1e-6);
nOfNodes = size(nodalXYZ, 2);

% Modify videoprefix, show number of nodes
videoprefix = strcat(videoprefix, "_Distance", num2str(Distance_To_Surface));
disp("Number of nodes in I1: ");
disp(sum(I1));

% Read nodal slip, slip rate and traction
Slip = h5read(faultFileName, '/vertex_fields/slip');
SlipRate = h5read(faultFileName, '/vertex_fields/slip_rate');
Traction = h5read(faultFileName, '/vertex_fields/traction');
Connection = h5read(faultFileName, '/topology/cells');

% Offset Traction by its initial values
if tractionOffsetFlag == 1
    Traction(1, :, :) = Traction(1, :, :) - Traction(1, :, 2) + 5.766727053863697e+06;
end

% Read nodal slip
% Input the first wire position
WirePos1 = [-0.025657; -0.014222; 0];
% VSstart = [0.020468, 0.011345, 0]';
% VSend = [0.079942, 0.044312, 0]';

VSstart = [0.006354, 0.003522, 0]';
VSend = [0.063204, 0.035034, 0]';
VSregion = 1e3 * [norm(VSstart - WirePos1, 2), norm(VSend - WirePos1, 2)];
% VSregion = [50, 120];

nOf2DNodes = sum(I1);

FaultX = zeros(1, nOf2DNodes);
surfaceNodesXYZ = nodalXYZ(:, I1);
surfaceTraction = Traction(1, I1, :);
surfaceNormal = Traction(3, I1, :);
surfaceTraction = reshape(surfaceTraction, size(surfaceTraction, 2), size(surfaceTraction, 3));
surfaceNormal = reshape(surfaceNormal, size(surfaceNormal, 2), size(surfaceNormal, 3));

[~, I] = sort(surfaceNodesXYZ(1,:));
surfaceTraction = surfaceTraction(I, :);
surfaceNormal = surfaceNormal(I, :);


% Fault Range
xrange = [-100, 150];
Vrange = [0, 20];

% Framerate for videos
framerate = 8;

for i = 1:1:nOf2DNodes
    FaultX(i) = norm(surfaceNodesXYZ(1:2, i) - WirePos1(1:2), 2) * sign(surfaceNodesXYZ(1, i) - WirePos1(1));
end
FaultX = FaultX(I);

% 2D node coordinates
nodalXYZ2D = zeros(2, nOfNodes);
for i = 1:1:nOfNodes
    nodalXYZ2D(1, i) = sign(nodalXYZ(1, i) - WirePos1(1)) ...
        * norm(nodalXYZ(1:2, i) - WirePos1(1:2), 2);
    nodalXYZ2D(2, i) = nodalXYZ(3, i);
end

% Magnitude of slip rate
slipRateMag = zeros(nOfNodes, nOfTimeSteps);
surfaceSlipRateMag = zeros(nOf2DNodes, nOfTimeSteps);
surfaceSlipRate = zeros(3, nOf2DNodes, nOfTimeSteps);

for t = 1:1:nOfTimeSteps
    surfaceSlipRate(:, :, t) = SlipRate(:, I1, t);
    for i = 1:1:nOf2DNodes
        surfaceSlipRateMag(i, t) = norm(surfaceSlipRate(:, i, t));
    end
    for i = 1:1:nOfNodes
        slipRateMag(i, t) = norm(SlipRate(:, i, t), 2);
    end
    surfaceSlipRateMag(:, t) = surfaceSlipRateMag(I, t);
end
figNo = 1;
%% Plot slip rate at different time
plotflag = true;
if plotflag == true
    dtstep = 10;
    tstep = 35;
    figure(figNo);
    labels = {};
    maxTimeStep = 65;
    while tstep <= maxTimeStep
        t = time(tstep);
        plot(1e3 * FaultX, surfaceSlipRateMag(:, tstep), 'linewidth', 2.0);
        xlim([0, 120]);
        ylim(Vrange);
        hold on; grid on;
        tstep = tstep + dtstep;
        labels{end + 1} = strcat('t = ', num2str(1e6 * t), ' \mu s ');
    end
    legend(labels, 'location', 'best');
    xlabel('Distance along the fault [mm]');
    ylabel('Slip Rate [m/s]');
    set(gca, 'FontSize', fontsize);
end
figNo = figNo + 1;

%% Save a X-T diagram plot of surface slip rate (normal range)
plotflag = true;
if plotflag == true
    fig = figure(figNo);
    Trange = [0, 115];
    Xrange = [-100, 150];
    fig.Position(3:4) = 1.5 * fig.Position(3:4);
    
    % Initialize names
    plotname = strcat(pwd, '/../plots/', videoprefix, '_X-TofSlipRate_normalRange.png');
    
    % Plot sliprate on X-T
    [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * FaultX);
    %subplot(2,2,iii)
    h = pcolor(Xsteps', Tsteps', surfaceSlipRateMag');
    shading interp;
    if VitoColorFlag == 1
        plotname = strcat(pwd, '/../Vitoplots/', videoprefix, '_X-TofSlipRate_normalRange.png');
        colormap(flipud(black_rainbow_shear_long));
    end
    hold on;
    xline(VSregion(1), 'r' ,'linewidth', 2.0);
    xline(VSregion(2), 'r' ,'linewidth', 2.0);
    text(VSregion(1)+ 5, 40, 'VS region', 'color', 'r', 'Fontsize', fontsize);
    % Add the wave speeds
    plot(cX, crY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, crY(2)+5, strcat('$c_r$ = 1.20 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 10, 'interpreter', 'latex');
    plot(cX, csY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, csY(2) - 1, strcat('$c_s$ = 1.28 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 10, 'interpreter', 'latex');
    plot(cX, cpY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, cpY(2), strcat('$c_p$ = 2.66 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 10, 'interpreter', 'latex');
    hold off;
    set(h, 'EdgeColor', 'None');
    c = colorbar;
    caxis([0, 12]);
    ylabel(c,'Slip rate [m/s]','FontName','Avenir','FontSize',fontsize);
    xlim(Xrange);
    ylim(Trange);
    xlabel('Distance along the fault [mm]');
    ylabel('Time [\mus]');
    title('X-T diagram of Slip rate');
    set(gca, 'FontSize', fontsize);
    
    % Save the figure
    print(figure(figNo) ,plotname, '-dpng', '-r500');
end
figNo = figNo + 1;

%% Save a X-T diagram plot of surface slip rate (Vito range)
plotflag = true;
if plotflag == true
    fig = figure(figNo);
    % Trange = [0, 150];
    Trange = [0, 115];
    Xrange = [-100, 150];
    fig.Position(3:4) = 1.5 * fig.Position(3:4);
    
    % Initialize names
    plotname = strcat(pwd, '/../plots/', videoprefix, '_X-TofSlipRate_VitoRange.png');
    
    % Plot sliprate on X-T
    [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * FaultX);
    %subplot(2,2,iii)
    h = pcolor(Xsteps', Tsteps', surfaceSlipRateMag');
    shading interp;
    if VitoColorFlag == 1
        plotname = strcat(pwd, '/../Vitoplots/', videoprefix, '_X-TofSlipRate_VitoRange.png');
        colormap(flipud(black_rainbow_shear_long));
    end
    hold on;
    xline(VSregion(1), 'r' ,'linewidth', 2.0);
    xline(VSregion(2), 'r' ,'linewidth', 2.0);
    text(VSregion(1)+ 5, 40, 'VS region', 'color', 'r', 'Fontsize', fontsize);
    % Add the wave speeds
    plot(cX, crY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, crY(2)+5, strcat('$c_r$ = 1.20 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 10, 'interpreter', 'latex');
    plot(cX, csY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, csY(2) - 1, strcat('$c_s$ = 1.28 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 10, 'interpreter', 'latex');
    plot(cX, cpY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, cpY(2), strcat('$c_p$ = 2.66 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 10, 'interpreter', 'latex');
    hold off;
    set(h, 'EdgeColor', 'None');
    c = colorbar;
    caxis([0, 2]);
    ylabel(c,'Slip rate [m/s]','FontName','Avenir','FontSize',fontsize);
    xlim(Xrange);
    ylim(Trange);
    xlabel('Distance along the fault [mm]');
    ylabel('Time [\mus]');
    title('X-T diagram of Slip rate');
    set(gca, 'FontSize', fontsize);
    
    % Save the figure
    print(figure(figNo) ,plotname, '-dpng', '-r500');
end
figNo = figNo + 1;

%% Save a X-T diagram plot of surface slip rate (normal range, log scale)
plotflag = true;
if plotflag == true
    fig = figure(figNo);
    % Trange = [0, 150];
    Trange = [0, 115];
    Xrange = [-100, 150];
    fig.Position(3:4) = 1.5 * fig.Position(3:4);
    
    % Initialize names
    plotname = strcat(pwd, '/../plots/', videoprefix, '_X-TofSlipRate_normalRangeLog.png');
    
    % Plot sliprate on X-T
    [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * FaultX);
    %subplot(2,2,iii)
    h = pcolor(Xsteps', Tsteps', log10(surfaceSlipRateMag)');
    shading interp;
    if VitoColorFlag == 1
        plotname = strcat(pwd, '/../Vitoplots/', videoprefix, '_X-TofSlipRate_normalRangeLog.png');
        colormap(flipud(black_rainbow_shear_long));
    end
    hold on;
    xline(VSregion(1), 'r' ,'linewidth', 2.0);
    xline(VSregion(2), 'r' ,'linewidth', 2.0);
    text(VSregion(1)+ 5, 40, 'VS region', 'color', 'r', 'Fontsize', fontsize);
    % Add the wave speeds
    plot(cX, crY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, crY(2)+5, strcat('$c_r$ = 1.20 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 10, 'interpreter', 'latex');
    plot(cX, csY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, csY(2) - 1, strcat('$c_s$ = 1.28 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 10, 'interpreter', 'latex');
    plot(cX, cpY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, cpY(2), strcat('$c_p$ = 2.66 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 10, 'interpreter', 'latex');
    hold off;
    set(h, 'EdgeColor', 'None');
    c = colorbar;
    % caxis([0, 12]);
    ylabel(c,'Log slip rate [m/s]','FontName','Avenir','FontSize',fontsize);
    xlim(Xrange);
    ylim(Trange);
    xlabel('Distance along the fault [mm]');
    ylabel('Time [\mus]');
    title('X-T diagram of Slip rate');
    set(gca, 'FontSize', fontsize);
    
    % Save the figure
    print(figure(figNo) ,plotname, '-dpng', '-r500');
end
figNo = figNo + 1;

%% Save a X-T diagram plot of surface slip rate (only observing window)
plotflag = true;
if plotflag == true
    fig = figure(figNo);
    % Trange = [0, 150];
    Trange = [30, 110];
    Xrange = [VSregion(1), VSregion(1) + 45];
    fig.Position(3:4) = 1.5 * fig.Position(3:4);
    
    % Initialize names
    plotname = strcat(pwd, '/../plots/', videoprefix, '_X-TofSlipRate_window.png');
    
    % Plot sliprate on X-T
    [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * FaultX);
    %subplot(2,2,iii)
    h = pcolor(Xsteps', Tsteps', surfaceSlipRateMag');
    shading interp;
    if VitoColorFlag == 1
        plotname = strcat(pwd, '/../Vitoplots/', videoprefix, '_X-TofSlipRate_window.png');
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
    caxis([0, 2]);
    ylabel(c,'Slip rate [m/s]','FontName','Avenir','FontSize',fontsize);
    xlim(Xrange);
    ylim(Trange);
    xlabel('Distance along the fault [mm]');
    ylabel('Time [\mus]');
    title('X-T diagram of Slip rate');
    set(gca, 'FontSize', fontsize);
    
    % Save the figure
    print(figure(figNo) ,plotname, '-dpng', '-r500');
end
figNo = figNo + 1;

%% Save a X-T diagram plot of surface slip rate (observing window, smaller range)
plotflag = true;
if plotflag == true
    fig = figure(figNo);
    % Trange = [0, 150];
    Trange = [30, 110];
    Xrange = [VSregion(1), VSregion(1) + 45];
    fig.Position(3:4) = 1.5 * fig.Position(3:4);
    
    % Initialize names
    plotname = strcat(pwd, '/../plots/', videoprefix, '_X-TofSlipRate_window_small.png');
    
    % Plot sliprate on X-T
    [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * FaultX);
    %subplot(2,2,iii)
    h = pcolor(Xsteps', Tsteps', surfaceSlipRateMag');
    shading interp;
    if VitoColorFlag == 1
        plotname = strcat(pwd, '/../Vitoplots/', videoprefix, '_X-TofSlipRate_window_small.png');
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
    title('X-T diagram of Slip rate');
    set(gca, 'FontSize', fontsize);
    
    % Save the figure
    print(figure(figNo) ,plotname, '-dpng', '-r500');
end
figNo = figNo + 1;

%% Save a X-T diagram plot of shear stress (normal range)
plotflag = true;
if plotflag == true
    fig = figure(figNo);
    Trange = [0, 115];
    Xrange = [-100, 150];
    fig.Position(3:4) = 1.5 * fig.Position(3:4);
    
    % Initialize names
    plotname = strcat(pwd, '/../plots/', videoprefix, '_X-TofShearStress_normalRange.png');
    
    % Plot sliprate on X-T
    [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * FaultX);
    %subplot(2,2,iii)
    h = pcolor(Xsteps', Tsteps', (surfaceTraction / 1.0e6)');
    shading interp;
    if VitoColorFlag == 1
        plotname = strcat(pwd, '/../Vitoplots/', videoprefix, '_X-TofShearStress_normalRange.png');
        colormap(black_rainbow_plus_long);
    end
    hold on;
    xline(VSregion(1), 'r' ,'linewidth', 2.0);
    xline(VSregion(2), 'r' ,'linewidth', 2.0);
    text(VSregion(1)+ 5, 40, 'VS region', 'color', 'r', 'Fontsize', fontsize);
    
    cX = [60, 80];
    crY = [10, (cX(2) - cX(1)) * 1e3 / cr + 10];
    csY = [10, (cX(2) - cX(1)) * 1e3 / cs + 10];
    cpY = [10, (cX(2) - cX(1)) * 1e3 / cp + 10];
    % Add the wave speeds
    plot(cX, crY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, crY(2)+5, strcat('$c_r$ = 1.20 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 10, 'interpreter', 'latex');
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

%% Save a X-T diagram plot of shear stress (only observing window)
plotflag = true;
if plotflag == true
    fig = figure(figNo);
    % Trange = [0, 150];
    Trange = [30, 110];
    Xrange = [VSregion(1), VSregion(1) + 45];
    fig.Position(3:4) = 1.5 * fig.Position(3:4);
    
    % Initialize names
    plotname = strcat(pwd, '/../plots/', videoprefix, '_X-TofShearstress_window.png');
    
    % Plot sliprate on X-T
    [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * FaultX);
    %subplot(2,2,iii)
    h = pcolor(Xsteps', Tsteps', (surfaceTraction ./ 1e6)');
    shading interp;
    if VitoColorFlag == 1
        plotname = strcat(pwd, '/../Vitoplots/', videoprefix, '_X-TofShearStress_window.png');
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

%% Save a X-T diagram plot of normal stress (only observing window)
plotflag = true;
if plotflag == true
    fig = figure(figNo);
    % Trange = [0, 150];
    Trange = [30, 110];
    Xrange = [VSregion(1), VSregion(1) + 45];
    fig.Position(3:4) = 1.5 * fig.Position(3:4);
    
    % Initialize names
    plotname = strcat(pwd, '/../plots/', videoprefix, '_X-TofNormalstress_window.png');
    
    % Plot sliprate on X-T
    [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * FaultX);
    %subplot(2,2,iii)
    h = pcolor(Xsteps', Tsteps', -(surfaceNormal ./ 1e6)');
    shading interp;
    if VitoColorFlag == 1
        plotname = strcat(pwd, '/../Vitoplots/', videoprefix, '_X-TofNormalStress_window.png');
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
    caxis([5, 15]);
    ylabel(c,'Shear stress [MPa]','FontName','Avenir','FontSize',fontsize);
    xlim(Xrange);
    ylim(Trange);
    xlabel('Distance along the fault [mm]');
    ylabel('Time [\mus]');
    title('X-T diagram of Normal stress');
    set(gca, 'FontSize', fontsize);
    
    % Save the figure
    print(figure(figNo) ,plotname, '-dpng', '-r500');
end
figNo = figNo + 1;