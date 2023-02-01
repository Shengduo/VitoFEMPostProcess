function writeVideosTet_function(videoprefix, pre_time)
    % Read results from hdf5 files.
    load('BRColorScale.mat');
    VitoColorFlag = true;
    tractionOffsetFlag = false;
    yToxRatio = 2;
    
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
    xrange = [-100, 150];
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
        nodalXYZ2D(2, i) = 5e-3 - nodalXYZ(3, i);
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
        fig = figure(figNo);
        fig.Position = [1000, 597, 2800/2, 1680/2];
        xrange = 1e3 * [min(FaultX), max(FaultX)];
        yrange = [0, 10];

        % Initialize names
        videoname = strcat(videoprefix, num2str(yToxRatio), '_sliprateMag.avi');

        % Initialize video
        myVideo = VideoWriter(strcat('../PaperVideos/', videoname), 'Motion JPEG AVI');
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
            xlabel('Distance along the fault [mm]', 'interpreter', 'latex');
            ylabel({'Depth into', 'surface [mm]'}, 'interpreter', 'latex');
            title(['Time = ', num2str(1e6 * time(i), '%.2f'), ' [$\mathrm{\mu s}$]'], 'interpreter', 'latex');
            set(gca, 'FontSize', fontsize);
            daspect([yToxRatio, 1, 1]);
            % Write the video
            set(gcf, 'color', 'w');
            % fig.Position = [1000, 597, 2800/2, 1680/2];
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
        fig.Position = [1000, 597, 2800/2, 1680/4];
        xrange = 1e3 * [min(FaultX), max(FaultX)];
        yrange = [0, 10];

        % Initialize names
        videoname = strcat(videoprefix, num2str(yToxRatio), '_shearTrac.avi');

        % Initialize video
        myVideo = VideoWriter(strcat('../PaperVideos/', videoname),'Motion JPEG AVI');
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
            colormap(black_rainbow_plus_long);
            hold on;
            xline(VSregion(1), 'r' ,'linewidth', 2.0);
            xline(VSregion(2), 'r' ,'linewidth', 2.0);
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
            xlabel('Distance along the fault [mm]', 'interpreter', 'latex');
            ylabel({'Depth into', 'surface [mm]'}, 'interpreter', 'latex');
            title(['Time = ', num2str(1e6 * time(i), '%.2f'), ' [$\mathrm{\mu s}$]'], 'interpreter', 'latex');
            set(gca, 'FontSize', fontsize);
            daspect([yToxRatio, 1, 1]);

            % Write the video
            set(gcf, 'color', 'w');
            fig.Position = [1000, 597, 2800/2, 1680/4];
            frame = getframe(gcf);
            writeVideo(myVideo, frame);
        end
        close(myVideo);
    end
    figNo = figNo + 1;

    %% Save a video of normal stress
    videoflag = false;
    if videoflag == true
        fig = figure(figNo);
        fig.Position = [1000, 597, 2800/2, 1680/2];
        xrange = 1e3 * [min(FaultX), max(FaultX)];
        yrange = [-5, 5];

        % Initialize names
        videoname = strcat(videoprefix, num2str(yToxRatio), '_normalTrac.avi');

        % Initialize video
        myVideo = VideoWriter(strcat('../PaperVideos/', videoname), 'Motion JPEG AVI');
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
            hold on;
            xline(VSregion(1), 'r' ,'linewidth', 2.0);
            xline(VSregion(2), 'r' ,'linewidth', 2.0);
            text(56 - 1e3 * norm(VSstart - WirePos1, 2), -2, 'Region 2', 'color', 'r', 'Fontsize', fontsize, 'interpreter', 'latex');
            hold off;
            c = colorbar;
            caxis([-12, 0]);
            c.Ticks = [-12, -6,  0];
            ylabel(c,{'Normal Stress', '[MPa]'},'FontName','Avenir','FontSize',fontsize, 'interpreter', 'latex');
            axis equal;
            grid on;
            xlim(xrange);
            ylim(yrange);
            xlabel('Distance along the fault [mm]', 'interpreter', 'latex');
            ylabel({'Interface', 'width [mm]'}, 'interpreter', 'latex');
            title(['Time = ', num2str(1e6 * time(i), '%.2f'), ' [$\mathrm{\mu s}$]'], 'interpreter', 'latex');
            set(gca, 'FontSize', fontsize);
            set(gca, 'ydir', 'reverse');
            daspect([yToxRatio, 1, 1]);

            % Write the video
            set(gcf, 'color', 'w');
            fig.Position = [1000, 597, 2800/2, 1680/2];
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
        myVideo = VideoWriter(strcat('../PaperVideos/', videoname), 'MPEG-4');
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
            fig.Position = [1000, 597, 2800/2, 1680/2];
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
        videoname = strcat(videoprefix, '_slip.avi');

        % Initialize video
        myVideo = VideoWriter(strcat('../PaperVideos/', videoname), 'MPEG-4');
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
            fig.Position = [1000, 597, 2800/2, 1680/2];
            frame = getframe(gcf);
            writeVideo(myVideo, frame);
        end
        close(myVideo);
    end
    figNo = figNo + 1;

    % %% Get the video of fault-parallel and fault-vertical velocity
    % % Read time
    % time = h5read(frontsurfFile, '/time');
    % time = reshape(time, [1, size(time, 3)]);
    % nOfTimeSteps = size(time, 2);
    % 
    % % Read node geometry
    % nodalXYZ = h5read(frontsurfFile, '/geometry/vertices');
    % nodalXYZ = nodalXYZ(1:2, :);
    % nOfNodes = size(nodalXYZ, 2);
    % 
    % % Read nodal slip, slip rate and traction
    % velocity = h5read(frontsurfFile, '/vertex_fields/velocity');
    % velocity = velocity(1:2, :, :);
    % Connection = h5read(frontsurfFile, '/topology/cells');
    % 
    % % Rotate into fault coordinate system
    % alpha = 29 / 180 * pi;
    % Q = [cos(alpha), sin(alpha); -sin(alpha), cos(alpha)];
    % 
    % % Rotate the velocity field
    % for i = 1:1:size(velocity, 3)
    %     velocity(:, :, i) = Q * velocity(:,:,i);
    % end
    % 
    % % Rotate nodal XYZ
    % nodalXYZ = Q * nodalXYZ;
    % 
    % % Set plot range
    % xrange = [floor(36.6 - norm(WirePos1)), ceil(36.6 - norm(WirePos1) + 47)];
    % yrange = [-15, 15];
    % crange = [-3, 3];
    % timeWindow = [25, 100];
    % 
    % %% Save a video of fault-parallel velocity field
    % videoflag = false;
    % if videoflag == true
    %     fig = figure(figNo);
    %     fig.Position = [1000, 597, 2240/4, 1680/4];
    %     
    %     % Initialize names
    %     videoname = strcat(videoprefix, '_prlV.mp4');
    %     
    %     % Initialize video
    %     myVideo = VideoWriter(strcat('../PaperVideos/', videoname), 'MPEG-4');
    %     myVideo.FrameRate = framerate;
    %     myVideo.Quality = 100;
    %     open(myVideo);
    %     
    %     % Shoot the video
    %     for i = 1:1:size(time, 2)
    %         if time(i) * 1e6 < timeWindow(1)
    %             continue;
    %         end
    %         p = patch(1e3 * nodalXYZ(1, :), 1e3 * nodalXYZ(2, :), ...
    %             velocity(1, :, i), 'EdgeColor', 'none', 'FaceColor', 'interp');
    %         p.Faces = Connection' + 1;
    %         yline(0, 'LineWidth', 2.0, 'color', 'w');
    %         c = colorbar;
    %         caxis(crange);
    %         ylabel(c,'Fault-parallel Velocity [m/s]','FontName','Avenir','FontSize',fontsize);
    %         axis equal;
    %         grid on;
    %         xlim(xrange);
    %         ylim(yrange);
    %         xlabel('X [mm]');
    %         ylabel('Y [mm]');
    %         title(strcat('Time = ', num2str(1e6 * time(i), '%2f'), ' [\mu s]'));
    %         set(gca, 'FontSize', fontsize);
    %         
    %         % Write the video
    %         frame = getframe(gcf);
    %         writeVideo(myVideo, frame);
    %         if time(i) * 1e6 > timeWindow(2)
    %             break;
    %         end
    %     end
    %     close(myVideo);
    % end
    % figNo = figNo + 1;
    % 
    % %% Save a video of fault-normal velocity field
    % videoflag = false;
    % if videoflag == true
    %     fig = figure(figNo);
    %     fig.Position = [1000, 597, 2240/4, 1680/4];
    %     
    %     % Initialize names
    %     videoname = strcat(videoprefix, '_nmlV.mp4');
    %     
    %     % Initialize video
    %     myVideo = VideoWriter(strcat('../PaperVideos/', videoname), 'MPEG-4');
    %     myVideo.FrameRate = framerate;
    %     myVideo.Quality = 100;
    %     open(myVideo);
    %     
    %     % Shoot the video
    %     for i = 1:1:size(time, 2)
    %         if time(i) * 1e6 < timeWindow(1)
    %             continue;
    %         end
    %         p = patch(1e3 * nodalXYZ(1, :), 1e3 * nodalXYZ(2, :), ...
    %             velocity(2, :, i), 'EdgeColor', 'none', 'FaceColor', 'interp');
    %         p.Faces = Connection' + 1;
    %         c = colorbar;
    %         caxis(crange);
    %         ylabel(c,'Fault-normal Velocity [m/s]','FontName','Avenir','FontSize',fontsize);
    %         axis equal;
    %         grid on;
    %         xlim(xrange);
    %         ylim(yrange);
    %         yline(0, 'LineWidth', 2.0, 'color', 'w');
    %         xlabel('X [mm]');
    %         ylabel('Y [mm]');
    %         title(strcat('Time = ', num2str(1e6 * time(i), '%2f'), ' [\mu s]'));
    %         set(gca, 'FontSize', fontsize);
    %         
    %         % Write the video
    %         frame = getframe(gcf);
    %         writeVideo(myVideo, frame);
    %         if time(i) * 1e6 > timeWindow(2)
    %             break;
    %         end
    %     end
    %     close(myVideo);
    % end
    % figNo = figNo + 1;
end