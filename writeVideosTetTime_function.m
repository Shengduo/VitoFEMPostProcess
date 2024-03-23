function writeVideosTet_function(videoprefix, xrange, pre_time, videoTimes, fmt, flag2D)
    % Read results from hdf5 files.
    load('BRColorScale.mat');
    VitoColorFlag = true;
    tractionOffsetFlag = false;
    yToxRatio = 2;
    
    X_total_range = [-121.599801068916, 107.070134864957];
    if xrange == 0.
        xrange = X_total_range;
    end
    xRatio = (xrange(2) - xrange(1)) / (X_total_range(2) - X_total_range(1));

    timeWindow = [pre_time * 1e6, 1000];
    % videoprefix = '1NPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.00224_-9_NULoad2dir0_duration120';
    % videoprefix = 'ViscoElastic_theta0.043';
    % videoprefix = 'DiffNULoadWithWallDRS1.5_8ModA0.016Load5_Vw2_fw0.1_theta0.036_8_NULoad2dir-1';
    faultFileName = strcat('../faultFiles/', videoprefix, '-fault.h5');
    h5disp(faultFileName);
    fontsize = 18;

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

    % Select the surface nodes
    I1 = (abs(nodalXYZ(3, :) - 0.005) <= 1e-7);
    nOfNodes = size(nodalXYZ, 2);

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
    VSregion = 1e3 * ([norm(VSstart - WirePos1, 2), norm(VSend - WirePos1, 2)] - norm(VSstart - WirePos1, 2));
    % VSregion = [50, 120];

    nOf2DNodes = sum(I1);

    FaultX = zeros(1, nOf2DNodes);
    surfaceNodesXYZ = nodalXYZ(:, I1);
    surfaceTraction = Traction(1, I1, :);
    surfaceTraction = reshape(surfaceTraction, size(surfaceTraction, 2), size(surfaceTraction, 3));
    [~, I] = sort(surfaceNodesXYZ(1,:));
    surfaceTraction = surfaceTraction(I, :);

    % Fault Range
%     xrange = [-100, 150];
    Vrange = [0, 20];

    % Framerate for videos
    framerate = 8;

    for i = 1:1:nOf2DNodes
        FaultX(i) = norm(surfaceNodesXYZ(1:2, i) - WirePos1(1:2), 2) * sign(surfaceNodesXYZ(1, i) - WirePos1(1));
    end
    FaultX = FaultX(I);
    FaultX = FaultX - norm(VSstart - WirePos1, 2);

    % 2D node coordinates
    nodalXYZ2D = zeros(2, nOfNodes);
    for i = 1:1:nOfNodes
        nodalXYZ2D(1, i) = sign(nodalXYZ(1, i) - VSstart(1)) ...
            * norm(nodalXYZ(1:2, i) - VSstart(1:2), 2);
        nodalXYZ2D(2, i) = nodalXYZ(3, i) - 5.e-3;
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

    %% Save a video of sliprate magnitude on the fault
    videoflag = false;
    if videoflag == true
        % Write the figures
        for idx = 1:1:size(videoTimes, 2)
            fig = figure(figNo);
            fig.Position = [1000, 597, 2400/2 * sqrt(xRatio), 250 / sqrt(xRatio)];

            % xrange = 1e3 * [min(FaultX), max(FaultX)];
            yrange = [-10, 0];

            [~, i] = min(abs(time - videoTimes(idx))); 

            p = patch(1e3 * nodalXYZ2D(1, :), 1e3 * nodalXYZ2D(2, :), ...
                slipRateMag(:, i), 'edgecolor', 'none');
            hold on;
            xline(VSregion(1), 'r' ,'linewidth', 2.0);
            xline(VSregion(2), 'r' ,'linewidth', 2.0);
            text(50, 4, strcat(num2str(1e6 * videoTimes(idx)), "$\ \mathrm{\mu s}$"), 'color', 'w', 'Fontsize', fontsize, 'interpreter', 'latex');
            text(56 - 1e3 * norm(VSstart - WirePos1, 2), -2, 'Region 2', 'color', 'r', 'Fontsize', fontsize, 'interpreter', 'latex');
            hold off;
            p.Faces = Connection' + 1;
            colormap(flipud(black_rainbow_shear_long));
            c = colorbar;
            caxis([0, 2]);
            c.Ticks = [0, 1, 2];
            ylabel(c,{'Slip rate', '[m/s]'},'FontName','Avenir','FontSize',fontsize,'interpreter', 'latex');
            axis equal;
            grid on;
            xlim(xrange);
            ylim(yrange);
            xlabel('$x_1\ [\mathrm{mm}]$', 'interpreter', 'latex');
            ylabel('$x_3\ [\mathrm{mm}]$', 'interpreter', 'latex');
            title(['Time = ', num2str(1e6 * time(i), '%.2f'), ' [$\mathrm{\mu s}$]'], 'interpreter', 'latex');
            set(gca, 'FontSize', fontsize);
            daspect([yToxRatio, 1, 1]);
            set(gca, 'YDir','reverse')

            % Write the video
            set(gcf, 'color', 'w');
            % fig.Position = [1000, 597, 2400/2, 1680/2];
            
            % Save the figure
            if fmt == "png"
                plotname = strcat(pwd, '/../PaperPlots_interp_snapshots/', videoprefix, '_V_', num2str(1e6 * videoTimes(idx)), '.png');
                print(figure(figNo), plotname, '-dpng', '-r250');
            else
                plotname = strcat(pwd, '/../PaperPlots_interp_snapshots/', videoprefix, '_V_', num2str(1e6 * videoTimes(idx)), '.eps');
                print('-vector', figure(figNo) ,plotname, '-depsc', '-r250');
            end
            close(fig);
        end
    end

    %% Save a video of shear stress
    videoflag = false;
    if videoflag == true
        % Write the figures
        for idx = 1:1:size(videoTimes, 2)
            fig = figure(figNo);
            
            fig.Position = [1000, 597, 2400/2 * sqrt(xRatio), 250 / sqrt(xRatio)];
            % xrange = 1e3 * [min(FaultX), max(FaultX)];
            yrange = [-10., 0.];
            [~, i] = min(abs(time - videoTimes(idx))); 
            p = patch(1e3 * nodalXYZ2D(1, :), 1e3 * nodalXYZ2D(2, :), ...
                Traction(1, :, i) / 1.0e6, 'edgecolor', 'none');
            p.Faces = Connection' + 1;
            colormap(black_rainbow_plus_long);
            hold on;
            xline(VSregion(1), 'r' ,'linewidth', 2.0);
            xline(VSregion(2), 'r' ,'linewidth', 2.0);
            text(50, 4, strcat(num2str(1e6 * videoTimes(idx)), "$\ \mathrm{\mu s}$"), 'color', 'w', 'Fontsize', fontsize, 'interpreter', 'latex');
            text(56 - 1e3 * norm(VSstart - WirePos1, 2), -2, 'Region 2', 'color', 'r', 'Fontsize', fontsize, 'interpreter', 'latex');
            hold off;
            c = colorbar;
            caxis([2, 10]);
            c.Ticks = [2, 6, 10];
            ylabel(c,{'Shear Stress', '[MPa]'},'FontName','Avenir','FontSize',fontsize, 'interpreter', 'latex');
            axis equal;
            grid on;
            xlim(xrange);
            ylim(yrange);
            xlabel('$x_1\ [\mathrm{mm}]$', 'interpreter', 'latex');
            ylabel('$x_3\ [\mathrm{mm}]$', 'interpreter', 'latex');
            title(['Time = ', num2str(1e6 * time(i), '%.2f'), ' [$\mathrm{\mu s}$]'], 'interpreter', 'latex');
            set(gca, 'FontSize', fontsize);
            daspect([yToxRatio, 1, 1]);
            set(gca, 'YDir','reverse')

            % Write the video
            set(gcf, 'color', 'w');
             % Save the figure
            if fmt == "png"
                plotname = strcat(pwd, '/../PaperPlots_interp_snapshots/', videoprefix, '_tau_', num2str(1e6 * videoTimes(idx)), '.png');
                print(figure(figNo), plotname, '-dpng', '-r250');
            else
                plotname = strcat(pwd, '/../PaperPlots_interp_snapshots/', videoprefix, '_tau_', num2str(1e6 * videoTimes(idx)), '.eps');
                print('-vector', figure(figNo) ,plotname, '-depsc', '-r250');
            end
            close(fig);

        end
    end

    %% Save a video of 2D slip rate
    if flag2D == 1
        load("../matFiles/2DVideo.mat");

        % Write the figures
        for idx = 1:1:size(videoTimes, 2)
            fig = figure(figNo);
            
            fig.Position = [1000, 597, 2400/2 * sqrt(xRatio), 250 / sqrt(xRatio)];
            % xrange = 1e3 * [min(FaultX), max(FaultX)];
            yrange = [-10., 0.];
            [~, i] = min(abs(TwoDTime - videoTimes(idx))); 
            slipRate2Dto3D = interp1(TwoDNodalXYZ, TwoDSlipRateMag(:, i), nodalXYZ2D(1, :));

            p = patch(1e3 * nodalXYZ2D(1, :), 1e3 * nodalXYZ2D(2, :), ...
                slipRate2Dto3D, 'edgecolor', 'none');
            hold on;
            xline(VSregion(1), 'r' ,'linewidth', 2.0);
            xline(VSregion(2), 'r' ,'linewidth', 2.0);
            text(50, 4, strcat(num2str(1e6 * videoTimes(idx)), "$\ \mathrm{\mu s}$"), 'color', 'w', 'Fontsize', fontsize, 'interpreter', 'latex');
            text(56 - 1e3 * norm(VSstart - WirePos1, 2), -2, 'Region 2', 'color', 'r', 'Fontsize', fontsize, 'interpreter', 'latex');
            hold off;
            p.Faces = Connection' + 1;
            colormap(flipud(black_rainbow_shear_long));
            c = colorbar;
            caxis([0, 2]);
            c.Ticks = [0, 1, 2];
            ylabel(c,{'Slip rate', '[m/s]'},'FontName','Avenir','FontSize',fontsize,'interpreter', 'latex');
            axis equal;
            grid on;
            xlim(xrange);
            ylim(yrange);
            xlabel('$x_1\ [\mathrm{mm}]$', 'interpreter', 'latex');
            ylabel('$x_3\ [\mathrm{mm}]$', 'interpreter', 'latex');
            title(['Time = ', num2str(1e6 * TwoDTime(i), '%.2f'), ' [$\mathrm{\mu s}$]'], 'interpreter', 'latex');
            set(gca, 'FontSize', fontsize);
            daspect([yToxRatio, 1, 1]);
            set(gca, 'YDir','reverse')

            % Write the video
            set(gcf, 'color', 'w');
            % fig.Position = [1000, 597, 2400/2, 1680/2];
            
            % Save the figure
            if fmt == "png"
                plotname = strcat(pwd, '/../PaperPlots_interp_snapshots/', '2DV_ini-7', '_V_', num2str(1e6 * videoTimes(idx)), '.png');
                print(figure(figNo), plotname, '-dpng', '-r250');
            else
                plotname = strcat(pwd, '/../PaperPlots_interp_snapshots/', '2DV_ini-7', '_V_', num2str(1e6 * videoTimes(idx)), '.eps');
                print('-vector', figure(figNo) ,plotname, '-depsc', '-r250');
            end
            close(fig);
        end
    end
end