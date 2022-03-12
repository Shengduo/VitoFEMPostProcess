% Read results from hdf5 files.
clc,clear;
close all;
load('BRColorScale.mat');
VitoColorFlag = true;

% frontsurfFile = 'FineFHLoad15DRS1-frontsurf.h5';
videoprefix = 'VWFirst_DRS1_8ModA0.02Load4.75_Vww0.1_1.5_4_fw0.35_0.1_theta0.09';
faultFileName = strcat(videoprefix, '-fault.h5');
h5disp(faultFileName);
fontsize = 25;

% Read time
time = h5read(faultFileName, '/time');
time = reshape(time, [1, size(time, 3)]);
time = time - 10e-6;
nOfTimeSteps = size(time, 2);

% Several wave speed to show on the plot
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

% Select the surface nodes
I1 = (abs(nodalXYZ(3, :) - 0.005) <= 1e-7);
nOfNodes = size(nodalXYZ, 2);

% Read nodal slip, slip rate and traction
Slip = h5read(faultFileName, '/vertex_fields/slip');
SlipRate = h5read(faultFileName, '/vertex_fields/slip_rate');
Traction = h5read(faultFileName, '/vertex_fields/traction');
Connection = h5read(faultFileName, '/topology/cells');

% Read nodal slip
% Input the first wire position
WirePos1 = [-0.025657; -0.014222; 0];
% VSstart = [0.020468, 0.011345, 0]';
% VSend = [0.079942, 0.044312, 0]';

VSstart = [0.006354, 0.003522, 0]';
VSend = [0.058832, 0.032610, 0]';
VSregion = 1e3 * [norm(VSstart - WirePos1, 2), norm(VSend - WirePos1, 2)];
% VSregion = [50, 120];

nOf2DNodes = sum(I1);

FaultX = zeros(1, nOf2DNodes);
surfaceNodesXYZ = nodalXYZ(:, I1);
surfaceTraction = Traction(1, I1, :);
surfaceTraction = reshape(surfaceTraction, size(surfaceTraction, 2), size(surfaceTraction, 3));
[~, I] = sort(surfaceNodesXYZ(1,:));
surfaceTraction = surfaceTraction(I, :);

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

%% Save a video of surface slip rate
videoflag = false;
if videoflag == true
    figure(figNo);
    yrange = [0, 20];
    
    % Initialize names
    videoname = strcat(videoprefix, '_surface.mp4');
    
    % Initialize video
    myVideo = VideoWriter(strcat('Videos/', videoname), 'MPEG-4');
    myVideo.FrameRate = framerate;
    myVideo.Quality = 100;
    open(myVideo);
    
    % Shoot the video
    for i = 1:1:size(time, 2)
        plot(1e3 * FaultX, surfaceSlipRateMag(:, i), 'linewidth', 2.0);
        grid on;
        xlim(xrange);
        ylim(yrange);
        xlabel('Distance along the fault [mm]');
        ylabel('Slip Rate [m/s]');
        title(strcat('Time = ', num2str(1e6 * time(i), '%2f'), ' [\mu s]'));
        set(gca, 'FontSize', fontsize);
        
        % Write the video
        frame = getframe(gcf);
        writeVideo(myVideo, frame);
    end
    close(myVideo);
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
    plotname = strcat(pwd, '/plots/', videoprefix, '_X-TofSlipRate_normalRange.png');
    
    % Plot sliprate on X-T
    [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * FaultX);
    %subplot(2,2,iii)
    h = pcolor(Xsteps', Tsteps', surfaceSlipRateMag');
    shading interp;
    if VitoColorFlag == 1
        plotname = strcat(pwd, '/Vitoplots/', videoprefix, '_X-TofSlipRate_normalRange.png');
        colormap(flipud(black_rainbow_shear_long));
    end
    hold on;
    xline(VSregion(1), 'r' ,'linewidth', 2.0);
    xline(VSregion(2), 'r' ,'linewidth', 2.0);
    text(VSregion(1)+ 5, 40, 'VS region', 'color', 'r', 'Fontsize', fontsize);
    % Add the wave speeds
    plot(cX, crY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, crY(2)+5, strcat('$c_r$ = 1.20 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
    plot(cX, csY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, csY(2) - 1, strcat('$c_s$ = 1.28 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
    plot(cX, cpY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, cpY(2), strcat('$c_p$ = 2.66 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
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
    plotname = strcat(pwd, '/plots/', videoprefix, '_X-TofSlipRate_VitoRange.png');
    
    % Plot sliprate on X-T
    [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * FaultX);
    %subplot(2,2,iii)
    h = pcolor(Xsteps', Tsteps', surfaceSlipRateMag');
    shading interp;
    if VitoColorFlag == 1
        plotname = strcat(pwd, '/Vitoplots/', videoprefix, '_X-TofSlipRate_VitoRange.png');
        colormap(flipud(black_rainbow_shear_long));
    end
    hold on;
    xline(VSregion(1), 'r' ,'linewidth', 2.0);
    xline(VSregion(2), 'r' ,'linewidth', 2.0);
    text(VSregion(1)+ 5, 40, 'VS region', 'color', 'r', 'Fontsize', fontsize);
    % Add the wave speeds
    plot(cX, crY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, crY(2)+5, strcat('$c_r$ = 1.20 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
    plot(cX, csY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, csY(2) - 1, strcat('$c_s$ = 1.28 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
    plot(cX, cpY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, cpY(2), strcat('$c_p$ = 2.66 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
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
    plotname = strcat(pwd, '/plots/', videoprefix, '_X-TofSlipRate_normalRangeLog.png');
    
    % Plot sliprate on X-T
    [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * FaultX);
    %subplot(2,2,iii)
    h = pcolor(Xsteps', Tsteps', log10(surfaceSlipRateMag)');
    shading interp;
    if VitoColorFlag == 1
        plotname = strcat(pwd, '/Vitoplots/', videoprefix, '_X-TofSlipRate_normalRangeLog.png');
        colormap(flipud(black_rainbow_shear_long));
    end
    hold on;
    xline(VSregion(1), 'r' ,'linewidth', 2.0);
    xline(VSregion(2), 'r' ,'linewidth', 2.0);
    text(VSregion(1)+ 5, 40, 'VS region', 'color', 'r', 'Fontsize', fontsize);
    % Add the wave speeds
    plot(cX, crY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, crY(2)+5, strcat('$c_r$ = 1.20 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
    plot(cX, csY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, csY(2) - 1, strcat('$c_s$ = 1.28 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
    plot(cX, cpY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, cpY(2), strcat('$c_p$ = 2.66 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
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
    plotname = strcat(pwd, '/plots/', videoprefix, '_X-TofSlipRate_window.png');
    
    % Plot sliprate on X-T
    [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * FaultX);
    %subplot(2,2,iii)
    h = pcolor(Xsteps', Tsteps', surfaceSlipRateMag');
    shading interp;
    if VitoColorFlag == 1
        plotname = strcat(pwd, '/Vitoplots/', videoprefix, '_X-TofSlipRate_window.png');
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
    text(cX(2) + 4, crY(2)+2, strcat('$c_r$ = 1.20 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
    plot(cX, csY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, csY(2) - 1, strcat('$c_s$ = 1.28 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
    plot(cX, cpY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, cpY(2), strcat('$c_p$ = 2.66 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
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
    plotname = strcat(pwd, '/plots/', videoprefix, '_X-TofSlipRate_window_small.png');
    
    % Plot sliprate on X-T
    [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * FaultX);
    %subplot(2,2,iii)
    h = pcolor(Xsteps', Tsteps', surfaceSlipRateMag');
    shading interp;
    if VitoColorFlag == 1
        plotname = strcat(pwd, '/Vitoplots/', videoprefix, '_X-TofSlipRate_window_small.png');
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
    text(cX(2) + 4, crY(2)+2, strcat('$c_r$ = 1.20 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
    plot(cX, csY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, csY(2) - 1, strcat('$c_s$ = 1.28 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
    plot(cX, cpY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, cpY(2), strcat('$c_p$ = 2.66 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
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
    plotname = strcat(pwd, '/plots/', videoprefix, '_X-TofShearStress_normalRange.png');
    
    % Plot sliprate on X-T
    [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * FaultX);
    %subplot(2,2,iii)
    h = pcolor(Xsteps', Tsteps', (surfaceTraction / 1.0e6)');
    shading interp;
    if VitoColorFlag == 1
        plotname = strcat(pwd, '/Vitoplots/', videoprefix, '_X-TofShearStress_normalRange.png');
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
    text(cX(2) + 4, crY(2)+5, strcat('$c_r$ = 1.20 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
    plot(cX, csY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, csY(2) - 1, strcat('$c_s$ = 1.28 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
    plot(cX, cpY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, cpY(2), strcat('$c_p$ = 2.66 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
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
    plotname = strcat(pwd, '/plots/', videoprefix, '_X-TofShearstress_window.png');
    
    % Plot sliprate on X-T
    [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * FaultX);
    %subplot(2,2,iii)
    h = pcolor(Xsteps', Tsteps', (surfaceTraction ./ 1e6)');
    shading interp;
    if VitoColorFlag == 1
        plotname = strcat(pwd, '/Vitoplots/', videoprefix, '_X-TofShearStress_window.png');
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
    text(cX(2) + 4, crY(2)+2, strcat('$c_r$ = 1.20 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
    plot(cX, csY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, csY(2) - 1, strcat('$c_s$ = 1.28 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
    plot(cX, cpY, 'w', 'linewidth', 2.0);
    text(cX(2) + 4, cpY(2), strcat('$c_p$ = 2.66 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
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

timeWindow = [0, 150];
%% Save a video of sliprate magnitude on the fault
videoflag = true;
if videoflag == true
    fig = figure(figNo);
    fig.Position = [1000, 597, 2800/2, 1680/2];
    % xrange = [-100, 150];
    yrange = [-6, 6];
    
    % Initialize names
    videoname = strcat(videoprefix, '_sliprateMag.mp4');
    
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
        p = patch(1e3 * nodalXYZ2D(1, :), 1e3 * nodalXYZ2D(2, :), ...
            slipRateMag(:, i), 'edgecolor', 'none');
        hold on;
        xline(VSregion(1), 'r' ,'linewidth', 2.0);
        xline(VSregion(2), 'r' ,'linewidth', 2.0);
        text(VSregion(1)+ 5, -8, 'VS region', 'color', 'r', 'Fontsize', fontsize);
        hold off;
        p.Faces = Connection' + 1;
        c = colorbar;
        caxis([0, 12]);
        ylabel(c,'Slip rate [m/s]','FontName','Avenir','FontSize',fontsize);
        axis equal;
        grid on;
        xlim(xrange);
        ylim(yrange);
        xlabel('Distance along the fault [mm]');
        ylabel({'Thickness along', 'the fault [mm]'});
        title(strcat('Time = ', num2str(1e6 * time(i), '%2f'), ' [\mus]'));
        set(gca, 'FontSize', fontsize);
        
        % Write the video
        frame = getframe(gcf);
        writeVideo(myVideo, frame);
    end
    close(myVideo);
end
figNo = figNo + 1;

%% Save a video of shear stress
videoflag = true;
if videoflag == true
    fig = figure(figNo);
    fig.Position = [1000, 597, 2800/2, 1680/2];
    % xrange = [-100, 150];
    yrange = [-6, 6];
    
    % Initialize names
    videoname = strcat(videoprefix, '_shearTrac.mp4');
    
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
        p = patch(1e3 * nodalXYZ2D(1, :), 1e3 * nodalXYZ2D(2, :), ...
            Traction(1, :, i) / 1.0e6, 'edgecolor', 'none');
        p.Faces = Connection' + 1;
        hold on;
        xline(VSregion(1), 'r' ,'linewidth', 2.0);
        xline(VSregion(2), 'r' ,'linewidth', 2.0);
        text(VSregion(1)+ 5, -8, 'VS region', 'color', 'r', 'Fontsize', fontsize);
        hold off;
        c = colorbar;
        caxis([0, 10]);
        ylabel(c,'Shear Stress [MPa]','FontName','Avenir','FontSize',fontsize);
        axis equal;
        grid on;
        xlim(xrange);
        ylim(yrange);
        xlabel('Distance along the fault [mm]');
        ylabel({'Thickness along', 'the fault [mm]'});
        title(strcat('Time = ', num2str(1e6 * time(i), '%2f'), ' [\mus]'));
        set(gca, 'FontSize', fontsize);
        
        % Write the video
        frame = getframe(gcf);
        writeVideo(myVideo, frame);
    end
    close(myVideo);
end
figNo = figNo + 1;

%% Save a video of normal stress
videoflag = true;
if videoflag == true
    fig = figure(figNo);
    fig.Position = [1000, 597, 2240/2, 1680/2];
    % xrange = [-100, 150];
    yrange = [-6, 6];
    
    % Initialize names
    videoname = strcat(videoprefix, '_normalTrac.mp4');
    
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
        p = patch(1e3 * nodalXYZ2D(1, :), 1e3 * nodalXYZ2D(2, :), ...
            Traction(3, :, i) / 1.0e6, 'EdgeColor', 'none');
        p.Faces = Connection' + 1;
        c = colorbar;
        caxis([-15, 5]);
        ylabel(c,'Normal Stress [MPa]','FontName','Avenir','FontSize',fontsize);
        axis equal;
        grid on;
        xlim(xrange);
        ylim(yrange);
        xlabel('Distance along the fault [mm]');
        ylabel({'Thickness along', 'the fault [mm]'});
        title(strcat('Time = ', num2str(1e6 * time(i), '%2f'), ' [\mu s]'));
        set(gca, 'FontSize', fontsize);
        
        % Write the video
        frame = getframe(gcf);
        writeVideo(myVideo, frame);
    end
    close(myVideo);
end
figNo = figNo + 1;

%% Save a video of fault open
videoflag = false;
if videoflag == true
    fig = figure(figNo);
    fig.Position = [1000, 597, 2240/2, 1680/2];
    yrange = [-6, 6];
    
    % Initialize names
    videoname = strcat(videoprefix, '_open.mp4');
    
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
        p = patch(1e3 * nodalXYZ2D(1, :), 1e3 * nodalXYZ2D(2, :), ...
            1e3 * Slip(2, :, i), 'EdgeColor', 'none');
        
        p.Faces = Connection' + 1;
        c = colorbar;
        caxis([-20, 20]);
        ylabel(c,'Fault Opening [mm]','FontName','Avenir','FontSize',fontsize);
        axis equal;
        grid on;
        xlim([-100, 150]);
        ylim(yrange);
        xlabel('Distance along the fault [mm]');
        ylabel({'Thickness along', 'the fault [mm]'});
        title(strcat('Time = ', num2str(1e6 * time(i), '%2f'), ' [\mu s]'));
        set(gca, 'FontSize', fontsize);
        
        % Write the video
        frame = getframe(gcf);
        writeVideo(myVideo, frame);
    end
    close(myVideo);
end
figNo = figNo + 1;

%% Save a video of fault slip
videoflag = false;
if videoflag == true
    figure(figNo);
    yrange = [0, 20];
    
    % Initialize names
    videoname = strcat(videoprefix, '_slip.mp4');
    
    % Initialize video
    myVideo = VideoWriter(strcat('Videos/', videoname), 'MPEG-4');
    myVideo.FrameRate = framerate;
    myVideo.Quality = 100;
    open(myVideo);
    
    % Shoot the video
    for i = 1:1:size(time, 2)
        plot(1e3 * FaultX, surfaceSlipRateMag(:, i), 'linewidth', 2.0);
        grid on;
        xlim(xrange);
        ylim(yrange);
        xlabel('Distance along the fault [mm]');
        ylabel('Slip Rate [m/s]');
        title(strcat('Time = ', num2str(1e6 * time(i), '%2f'), ' [\mu s]'));
        set(gca, 'FontSize', fontsize);
        
        % Write the video
        frame = getframe(gcf);
        writeVideo(myVideo, frame);
    end
    close(myVideo);
end
figNo = figNo + 1;

%% Get the video of fault-parallel and fault-vertical velocity
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
Connection = h5read(frontsurfFile, '/topology/cells');

% Rotate into fault coordinate system
alpha = 29 / 180 * pi;
Q = [cos(alpha), sin(alpha); -sin(alpha), cos(alpha)];

% Rotate the velocity field
for i = 1:1:size(velocity, 3)
    velocity(:, :, i) = Q * velocity(:,:,i);
end

% Rotate nodal XYZ
nodalXYZ = Q * nodalXYZ;

% Set plot range
xrange = [floor(36.6 - norm(WirePos1)), ceil(36.6 - norm(WirePos1) + 47)];
yrange = [-15, 15];
crange = [-3, 3];
timeWindow = [25, 100];

%% Save a video of fault-parallel velocity field
videoflag = false;
if videoflag == true
    fig = figure(figNo);
    fig.Position = [1000, 597, 2240/4, 1680/4];
    
    % Initialize names
    videoname = strcat(videoprefix, '_prlV.mp4');
    
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
        p = patch(1e3 * nodalXYZ(1, :), 1e3 * nodalXYZ(2, :), ...
            velocity(1, :, i), 'EdgeColor', 'none', 'FaceColor', 'interp');
        p.Faces = Connection' + 1;
        yline(0, 'LineWidth', 2.0, 'color', 'w');
        c = colorbar;
        caxis(crange);
        ylabel(c,'Fault-parallel Velocity [m/s]','FontName','Avenir','FontSize',fontsize);
        axis equal;
        grid on;
        xlim(xrange);
        ylim(yrange);
        xlabel('X [mm]');
        ylabel('Y [mm]');
        title(strcat('Time = ', num2str(1e6 * time(i), '%2f'), ' [\mu s]'));
        set(gca, 'FontSize', fontsize);
        
        % Write the video
        frame = getframe(gcf);
        writeVideo(myVideo, frame);
        if time(i) * 1e6 > timeWindow(2)
            break;
        end
    end
    close(myVideo);
end
figNo = figNo + 1;

%% Save a video of fault-normal velocity field
videoflag = false;
if videoflag == true
    fig = figure(figNo);
    fig.Position = [1000, 597, 2240/4, 1680/4];
    
    % Initialize names
    videoname = strcat(videoprefix, '_nmlV.mp4');
    
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
        p = patch(1e3 * nodalXYZ(1, :), 1e3 * nodalXYZ(2, :), ...
            velocity(2, :, i), 'EdgeColor', 'none', 'FaceColor', 'interp');
        p.Faces = Connection' + 1;
        c = colorbar;
        caxis(crange);
        ylabel(c,'Fault-normal Velocity [m/s]','FontName','Avenir','FontSize',fontsize);
        axis equal;
        grid on;
        xlim(xrange);
        ylim(yrange);
        yline(0, 'LineWidth', 2.0, 'color', 'w');
        xlabel('X [mm]');
        ylabel('Y [mm]');
        title(strcat('Time = ', num2str(1e6 * time(i), '%2f'), ' [\mu s]'));
        set(gca, 'FontSize', fontsize);
        
        % Write the video
        frame = getframe(gcf);
        writeVideo(myVideo, frame);
        if time(i) * 1e6 > timeWindow(2)
            break;
        end
    end
    close(myVideo);
end
figNo = figNo + 1;