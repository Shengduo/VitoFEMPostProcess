function [X, Y, slipRates] = SlipRateAtDistInTheFault_function(videoprefix, Distance_To_Surface, pre_time, fmt)
    % This function process the faultfile and get the slip at certain location
    load('BRColorScale.mat');
    VitoColorFlag = true;

    faultFileName = strcat('../faultFiles/', videoprefix, '-fault.h5');
    % h5disp(faultFileName);
    fontsize = 25;
    yToxRatio = 30 / 37.5;
    ftBig = 10;
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

    % Read node geometry
    nodalXYZ = h5read(faultFileName, '/geometry/vertices');

    %% Get nodal information
    nOfNodes = size(nodalXYZ, 2);

    % Modify videoprefix, show number of nodes
    videoprefix = strcat(videoprefix);

    % Read nodal slip, slip rate and traction
    SlipRate = h5read(faultFileName, '/vertex_fields/slip_rate');
    connection = h5read(faultFileName, '/topology/cells');
    connection = connection + 1;
    traction = h5read(faultFileName, '/vertex_fields/traction');
    
    % Read nodal slip
    % Input the first wire position
    WirePos1 = [-0.025657; -0.014222; 0];
    FaultStart = -[0.100000, 0.055430, 0]';
    FaultEnd = [0.100000, 0.055430, 0]';
    VSstart = [0.006354, 0.003522, 0]';
    VSend = [0.063204, 0.035034, 0]';
    VSregion = 1e3 * ([norm(VSstart - WirePos1, 2), norm(VSend - WirePos1, 2)] - norm(VSstart - WirePos1, 2));
    % VSregion = [50, 120];

    % Calculate NodalXYZ2D 
    NodalXYZ2D = zeros(2, size(nodalXYZ, 2));
    for i = 1:1:size(nodalXYZ, 2)
        NodalXYZ2D(1, i) = sign(nodalXYZ(1, i) - VSstart(1)) * norm(nodalXYZ(1:2, i) - VSstart(1:2), 2);
        NodalXYZ2D(2, i) = nodalXYZ(3, i);
    end

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
%     VS_start_x = norm(VSstart - VSstart, 2);
%     VS_end_x = norm(VSend - VSstart, 2);
    Fault_start_x = -norm(FaultStart - VSstart, 2);
    Fault_end_x = norm(FaultEnd - VSstart, 2);
    
    % Camera specific values
    Pix_10mm = 81;  %Number of pixels in a grid of 10 mm ;  144 pixels in 2 mm -> 72 pixels in 1 mm
    pxsize=10/(Pix_10mm)*1e-3; % pixel size in m/pixel  (Pix_10mm pixels  in 10 mm)
    stepsize = 1;

    % First specific the x and y's
    % x_up = VS_start_x : stepsize * pxsize : VS_end_x;
    x_up = Fault_start_x : stepsize * pxsize : Fault_end_x;
    
    % Get the querying positions
    [Xq, Yq] = meshgrid(x_up, Distance_To_Surface);
    P = [Xq', Yq'];

    elementID = pointLocation(TR, P);
    SlipRateAtDist = zeros(size(Xq, 2), size(time, 2));
    ShearStressAtDist = zeros(size(Xq, 2), size(time, 2));
    for t = 1:1:size(time, 2)
        for ele = 1:1:size(elementID, 1)
            vel_x = scatteredInterpolant(NodalXYZ2D(:, connection(:, elementID(ele)))', squeeze(slipRateMag(connection(:, elementID(ele)), t)), 'natural');
            stress_x = scatteredInterpolant(NodalXYZ2D(:, connection(:, elementID(ele)))', squeeze(traction(1, connection(:, elementID(ele)), t))', 'natural');
            SlipRateAtDist(ele, t) = vel_x(x_up(ele), Distance_To_Surface);
            ShearStressAtDist(ele, t) = stress_x(x_up(ele), Distance_To_Surface);
        end
    end


    %% Save a X-T diagram plot of surface slip rate at stress_distance, z_location (observing window, smaller range)
    plotflag = true;
    if plotflag == true
        fig = figure(figNo);
        % Trange = [0, 150];
        Trange = [30, 110] + pre_time * 1e6;
        Xrange = [VSregion(1), VSregion(1) + 45];
        fig.Position(3:3) = 2.4 / 1.5 * fig.Position(3:3);
        fig.Position(4:4) = 2.4 * fig.Position(4:4);

        % Initialize names
        plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofSlipRate_window_small_surface_', num2str(Distance_To_Surface), '.eps');

        % Plot sliprate on X-T
        [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * x_up);
        h = pcolor(Xsteps', Tsteps', SlipRateAtDist');
        shading interp;

        if VitoColorFlag == 1
            plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofSlipRate_window_surface_', num2str(Distance_To_Surface), '.eps');
            colormap(flipud(black_rainbow_shear_long));
        end
        hold on;
        % xline(VSregion(1), 'r' ,'linewidth', 2.0);
        % xline(VSregion(2), 'r' ,'linewidth', 2.0);
        % text(VSregion(1)+ 5, 40, 'VS region', 'color', 'r', 'Fontsize', fontsize, 'interpreter', 'latex');
        % Add the wave speeds
        
        % Add x_3
        text(5, 105, strcat("$x_3$ = ", string(1000 * Distance_To_Surface - 5), "$\ \mathrm{mm}$"), 'color', 'w', 'Fontsize', fontsize + ftBig);
        
        cX = [50, 60] - 1e3 * norm(VSstart - WirePos1);
        crY = [60, (cX(2) - cX(1)) * 1e3 / cr + 60] + pre_time * 1e6;
        csY = [60, (cX(2) - cX(1)) * 1e3 / cs + 60] + pre_time * 1e6;
        cpY = [60, (cX(2) - cX(1)) * 1e3 / cp + 60] + pre_time * 1e6;

        plot(cX, crY, 'w', 'linewidth', 2.0);
        text(cX(2) + 1, crY(2) + 2, strcat('$c_r$ = 1.20 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 0, 'interpreter', 'latex');
        plot(cX, csY, 'w', 'linewidth', 2.0);
        text(cX(2) + 1, csY(2) - 1, strcat('$c_s$ = 1.28 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 0, 'interpreter', 'latex');
        plot(cX, cpY, 'w', 'linewidth', 2.0);
        text(cX(2) + 1, cpY(2) - 2, strcat('$c_p$ = 2.66 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 0, 'interpreter', 'latex');
        hold off;
        set(h, 'EdgeColor', 'None');
        c = colorbar;
        clim([0, 2]);
        ylabel(c,'Slip rate [m/s]','FontName','Avenir','FontSize',fontsize + ftBig, 'interpreter', 'latex');
        xlim(Xrange);
        ylim(Trange);
        xlabel('$x_1\ \mathrm{[mm]}$', 'interpreter', 'latex');
        ylabel('Time [$\mathrm{\mu s}$]', 'interpreter', 'latex');
        xticks(0:10:40);
        yticks(30:10:110);
        % title(['X-T diagram of Slip rate at Depth = ', ' ', num2str(5 - 1e3 * Distance_To_Surface, '%.0f'), ' ', '[mm]'], 'interpreter', 'latex');
        set(gca, 'FontSize', fontsize + ftBig);
        daspect([yToxRatio, 1, 1]);
        
        % Save the figure
        if fmt == 'png'
            plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofSlipRate_window_surface_', num2str(Distance_To_Surface), '.png');
            print(figure(figNo), plotname, '-dpng', '-r500');
        else
            print('-vector', figure(figNo) ,plotname, '-depsc', '-r500');
        end
    end
    figNo = figNo + 1;
    
    
    %% Save a X-T diagram plot of shear stress (only observing window)
    plotflag = true;
    if plotflag == true
        fig = figure(figNo);
        % Trange = [0, 150];
        Trange = [30, 110] + pre_time * 1e6;
        Xrange = [VSregion(1), VSregion(1) + 45];
        fig.Position(3:3) = 2.4 / 1.5 * fig.Position(3:3);
        fig.Position(4:4) = 2.4 * fig.Position(4:4);

        % Initialize names
        plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofShearStress_window_surface_', num2str(Distance_To_Surface), '.eps');

        % Plot sliprate on X-T
        [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * x_up);
        h = pcolor(Xsteps', Tsteps', (ShearStressAtDist ./ 1e6)');
        shading interp;
        
        if VitoColorFlag == 1
            plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofShearStress_window_surface_', num2str(Distance_To_Surface), '.eps');
            colormap(black_rainbow_plus_long);
        end
        hold on;
        % xline(VSregion(1), 'r' ,'linewidth', 2.0);
        % xline(VSregion(2), 'r' ,'linewidth', 2.0);
        % text(VSregion(1)+ 5, 40, 'VS region', 'color', 'r', 'Fontsize', fontsize);
        % Add the wave speeds

        cX = [50, 60] - 1e3 * norm(VSstart - WirePos1);
        crY = [40, (cX(2) - cX(1)) * 1e3 / cr + 40] + pre_time * 1e6;
        csY = [40, (cX(2) - cX(1)) * 1e3 / cs + 40] + pre_time * 1e6;
        cpY = [40, (cX(2) - cX(1)) * 1e3 / cp + 40] + pre_time * 1e6;

        plot(cX, crY, 'w', 'linewidth', 2.0);
        text(cX(2) + 1, crY(2) + 4, strcat('$c_r$ = 1.20 [km/s]'), 'color', 'w', 'Fontsize', fontsize + 0, 'interpreter', 'latex');
        plot(cX, csY, 'w', 'linewidth', 2.0);
        text(cX(2) + 1, csY(2) - 1, strcat('$c_s$ = 1.28 [km/s]'), 'color', 'w', 'Fontsize', fontsize + 0, 'interpreter', 'latex');
        plot(cX, cpY, 'w', 'linewidth', 2.0);
        text(cX(2) + 1, cpY(2) - 2, strcat('$c_p$ = 2.66 [km/s]'), 'color', 'w', 'Fontsize', fontsize + 0, 'interpreter', 'latex');
        hold off;
        
        % Add x_3
        text(5, 105, strcat("$x_3$ = ", string(1000 * Distance_To_Surface - 5), "$\ \mathrm{mm}$"), 'color', 'w', 'Fontsize', fontsize + ftBig);
        
        set(h, 'EdgeColor', 'None');
        c = colorbar;
        clim([2, 10]);
        ylabel(c,'Shear stress [MPa]','FontName','Avenir','FontSize',fontsize + ftBig, 'interpreter', 'latex');
        xlim(Xrange);
        ylim(Trange);
        xticks(0:10:40);
        yticks(30:10:110);
        xlabel('$x_1\ \mathrm{[mm]}$', 'interpreter', 'latex');
        ylabel('Time [$\mathrm{\mu s}$]', 'interpreter', 'latex');
        % title(['X-T of Shear Stress at  = Depth', ' ', num2str(5 - 1e3 * Distance_To_Surface, '%.0f'), ' ', '[mm]'], 'interpreter', 'latex');
        set(gca, 'FontSize', fontsize + ftBig);
        daspect([yToxRatio, 1, 1]);
        
        % Save the figure
        if fmt == 'png'
            plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofShearStress_window_surface_', num2str(Distance_To_Surface), '.png')
            print(figure(figNo), plotname, '-dpng', '-r500');
        else
            print('-vector', figure(figNo) ,plotname, '-depsc', '-r500');
        end
    end
    figNo = figNo + 1;
    
    
    %% Save a X-T diagram plot of surface slip rate at stress_distance, z_location (whole field, vito range)
    plotflag = true;
    if plotflag == true
        fig = figure(figNo);
        % Trange = [0, 150];
        Trange = [0, 110] + pre_time * 1e6;
        Xrange = 1e3 * [Fault_start_x, Fault_end_x];
        fig.Position(3:4) = 1.5 * fig.Position(3:4);

        % Initialize names
        plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofSlipRate_whole_surface_', num2str(Distance_To_Surface), '.eps');

        % Plot sliprate on X-T
        [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * x_up);
        h = pcolor(Xsteps', Tsteps', SlipRateAtDist');
        shading interp;

        if VitoColorFlag == 1
            plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofSlipRate_whole_Vito_surface_', num2str(Distance_To_Surface), '.eps');
            colormap(flipud(black_rainbow_shear_long));
        end
        hold on;
        xline(VSregion(1), 'r' ,'linewidth', 2.0);
        xline(VSregion(2), 'r' ,'linewidth', 2.0);
        text(VSregion(1)+ 8, 38, 'Region 2', 'color', 'r', 'Fontsize', fontsize);
        
        % Add x_3
        text(-110, 100, strcat("$x_3$ = ", string(1000 * Distance_To_Surface - 5), "$\ \mathrm{mm}$"), 'color', 'w', 'Fontsize', fontsize);
        % Add the wave speeds
        
        % Several wave speed to show on the aplot
        cp = 2662.4;
        cs = 1279;
        nu = 0.35;
        cr = (0.874 + 0.196 * nu - 0.043 * nu^2 - 0.055 * nu^3) * cs;
        cX = [60, 80] - 15;
        cX = cX - norm(VSstart - WirePos1, 2) * 1e3;
        crY = [10, (cX(2) - cX(1)) * 1e3 / cr + 10];
        csY = [10, (cX(2) - cX(1)) * 1e3 / cs + 10];
        cpY = [10, (cX(2) - cX(1)) * 1e3 / cp + 10];
        
        % Add the wave speeds
        plot(cX, crY, 'w', 'linewidth', 2.0);
        text(cX(2) + 1, crY(2)+5, strcat('$c_r$ = 1.20 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
        plot(cX, csY, 'w', 'linewidth', 2.0);
        text(cX(2) + 1, crY(2) - 1, strcat('$c_s$ = 1.28 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
        plot(cX, cpY, 'w', 'linewidth', 2.0);
        text(cX(2) + 1, cpY(2) - 2, strcat('$c_p$ = 2.66 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
        hold off;
        
        set(h, 'EdgeColor', 'None');
        c = colorbar;
        c2 = 15;
        clim([0, c2]);
        ylabel(c,'Slip rate [m/s]','FontName','Avenir','FontSize',fontsize, 'interpreter', 'latex');
        xlim(Xrange);
        ylim(Trange);
        title('X-T diagram of slip rate')
        xlabel('$x_1\ \mathrm{[mm]}$', 'interpreter', 'latex');
        ylabel('Time [$\mathrm{\mu s}$]', 'interpreter', 'latex');
        % title(['X-T diagram of Slip rate at Depth = ', ' ', num2str(5 - 1e3 * Distance_To_Surface, '%.0f'), ' ', '[mm]'], 'interpreter', 'latex');
        set(gca, 'FontSize', fontsize);
        % daspect([yToxRatio, 1, 1]);
        if c2 == 15
            plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofSlipRate_whole_Vito_surface_', num2str(Distance_To_Surface), 'full.eps');
        end
        % Save the figure
        if fmt == 'png'
            plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofSlipRate_whole_Vito_surface_', num2str(Distance_To_Surface), 'full.png');
            print(figure(figNo), plotname, '-dpng', '-r500');
        else
            print('-vector', figure(figNo) ,plotname, '-depsc', '-r500');
        end
    end
    figNo = figNo + 1;
    
    %% Save a X-T diagram plot of surface slip rate at stress_distance, z_location (whole field, vito range)
    plotflag = true;
    if plotflag == true
        fig = figure(figNo);
        % Trange = [0, 150];
        Trange = [0, 110] + pre_time * 1e6;
        Xrange = 1e3 * [Fault_start_x, Fault_end_x];
        fig.Position(3:4) = 1.5 * fig.Position(3:4);

        % Initialize names
        plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofSlipRate_whole_surface_', num2str(Distance_To_Surface), '.eps');

        % Plot sliprate on X-T
        [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * x_up);
        h = pcolor(Xsteps', Tsteps', SlipRateAtDist');
        shading interp;

        if VitoColorFlag == 1
            plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofSlipRate_whole_Vito_surface_', num2str(Distance_To_Surface), '.eps');
            colormap(flipud(black_rainbow_shear_long));
        end
        hold on;
        xline(VSregion(1), 'r' ,'linewidth', 2.0);
        xline(VSregion(2), 'r' ,'linewidth', 2.0);
        text(VSregion(1)+ 8, 38, 'Region 2', 'color', 'r', 'Fontsize', fontsize);
        
        % Add x_3
        text(-110, 100, strcat("$x_3$ = ", string(1000 * Distance_To_Surface - 5), "$\ \mathrm{mm}$"), 'color', 'w', 'Fontsize', fontsize);
        % Add the wave speeds
        
        % Several wave speed to show on the aplot
        cp = 2662.4;
        cs = 1279;
        nu = 0.35;
        cr = (0.874 + 0.196 * nu - 0.043 * nu^2 - 0.055 * nu^3) * cs;
        cX = [60, 80] - 15;
        cX = cX - norm(VSstart - WirePos1, 2) * 1e3;
        crY = [10, (cX(2) - cX(1)) * 1e3 / cr + 10];
        csY = [10, (cX(2) - cX(1)) * 1e3 / cs + 10];
        cpY = [10, (cX(2) - cX(1)) * 1e3 / cp + 10];
        
        % Add the wave speeds
        plot(cX, crY, 'w', 'linewidth', 2.0);
        text(cX(2) + 1, crY(2)+5, strcat('$c_r$ = 1.20 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
        plot(cX, csY, 'w', 'linewidth', 2.0);
        text(cX(2) + 1, crY(2) - 1, strcat('$c_s$ = 1.28 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
        plot(cX, cpY, 'w', 'linewidth', 2.0);
        text(cX(2) + 1, cpY(2) - 2, strcat('$c_p$ = 2.66 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
        hold off;
        
        set(h, 'EdgeColor', 'None');
        c = colorbar;
        c2 = 2;
        clim([0, c2]);
        ylabel(c,'Slip rate [m/s]','FontName','Avenir','FontSize',fontsize, 'interpreter', 'latex');
        xlim(Xrange);
        ylim(Trange);
        title('X-T diagram of slip rate')
        xlabel('$x_1\ \mathrm{[mm]}$', 'interpreter', 'latex');
        ylabel('Time [$\mathrm{\mu s}$]', 'interpreter', 'latex');
        % title(['X-T diagram of Slip rate at Depth = ', ' ', num2str(5 - 1e3 * Distance_To_Surface, '%.0f'), ' ', '[mm]'], 'interpreter', 'latex');
        set(gca, 'FontSize', fontsize);
        % daspect([yToxRatio, 1, 1]);
        if c2 == 15
            plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofSlipRate_whole_Vito_surface_', num2str(Distance_To_Surface), 'full.eps');
        end
        % Save the figure
        if fmt == 'png'
            plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofSlipRate_whole_Vito_surface_', num2str(Distance_To_Surface), '.png');
            print(figure(figNo), plotname, '-dpng', '-r500');
        else
            print('-vector', figure(figNo) ,plotname, '-depsc', '-r500');
        end
    end
    figNo = figNo + 1;
    
    % Store the slip rates data
    X = Xsteps'; Y = Tsteps'; slipRates = SlipRateAtDist';
    
    %% Save a X-T diagram plot of shear stress (whole field)
    plotflag = true;
    if plotflag == true
        fig = figure(figNo);
        % Trange = [0, 150];
        Trange = [0, 110] + pre_time * 1e6;
        Xrange = 1e3 * [Fault_start_x, Fault_end_x];
        fig.Position(3:4) = 1.5 * fig.Position(3:4);

        % Initialize names
        plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofShearStress_whole_surface_', num2str(Distance_To_Surface), '.eps');

        % Plot sliprate on X-T
        [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * x_up);
        h = pcolor(Xsteps', Tsteps', (ShearStressAtDist ./ 1e6)');
        shading interp;
        
        if VitoColorFlag == 1
            plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofShearStress_whole_Vito_surface_', num2str(Distance_To_Surface), '.eps');
            colormap(black_rainbow_plus_long);
        end
        hold on;
        xline(VSregion(1), 'r' ,'linewidth', 2.0);
        xline(VSregion(2), 'r' ,'linewidth', 2.0);
        text(VSregion(1)+ 8, 38, 'Region 2', 'color', 'r', 'Fontsize', fontsize);
        
        % Add x_3
        text(-110, 100, strcat("$x_3$ = ", string(1000 * Distance_To_Surface - 5), "$\ \mathrm{mm}$"), 'color', 'w', 'Fontsize', fontsize);
        
        % Add the wave speeds
        cp = 2662.4;
        cs = 1279;
        nu = 0.35;
        cr = (0.874 + 0.196 * nu - 0.043 * nu^2 - 0.055 * nu^3) * cs;
        cX = [60, 80] - 15;
        cX = cX - norm(VSstart - WirePos1, 2) * 1e3;
        crY = [10, (cX(2) - cX(1)) * 1e3 / cr + 10];
        csY = [10, (cX(2) - cX(1)) * 1e3 / cs + 10];
        cpY = [10, (cX(2) - cX(1)) * 1e3 / cp + 10];
        
        % Add the wave speeds
        plot(cX, crY, 'w', 'linewidth', 2.0);
        text(cX(2) + 1, crY(2)+5, strcat('$c_r$ = 1.20 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
        plot(cX, csY, 'w', 'linewidth', 2.0);
        text(cX(2) + 1, crY(2) - 1, strcat('$c_s$ = 1.28 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
        plot(cX, cpY, 'w', 'linewidth', 2.0);
        text(cX(2) + 1, cpY(2) - 2, strcat('$c_p$ = 2.66 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
        hold off;
        
        set(h, 'EdgeColor', 'None');
        c = colorbar;
        clim([2, 10]);
        ylabel(c,'Shear stress [MPa]','FontName','Avenir','FontSize',fontsize, 'interpreter', 'latex');
        xlim(Xrange);
        ylim(Trange);
        xlabel('$x_1\ \mathrm{[mm]}$', 'interpreter', 'latex');
        ylabel('Time [$\mathrm{\mu s}$]', 'interpreter', 'latex');
        title('X-T diagram of shear stress');
        % title(['X-T of Shear Stress at  = Depth', ' ', num2str(5 - 1e3 * Distance_To_Surface, '%.0f'), ' ', '[mm]'], 'interpreter', 'latex');
        set(gca, 'FontSize', fontsize);
        % daspect([yToxRatio, 1, 1]);
        
        % Save the figure
        if fmt == 'png'
            plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofShearStress_whole_Vito_surface_', num2str(Distance_To_Surface), '.png');
            print(figure(figNo), plotname, '-dpng', '-r500');
        else
            print('-vector', figure(figNo) ,plotname, '-depsc', '-r500');
        end
    end
    figNo = figNo + 1;

    %% Save a X-T diagram plot of surface slip rate at stress_distance, z_location (whole field, log scale)
    plotflag = true;
    if (plotflag == true) && (Distance_To_Surface == 0.) 
        fig = figure(figNo);
        % Trange = [0, 150];
        Trange = [0, 110] + pre_time * 1e6;
        Xrange = 1e3 * [Fault_start_x, Fault_end_x];
        fig.Position(3:4) = 1.5 * fig.Position(3:4);

        % Initialize names
        plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofLogSlipRate_whole_surface', num2str(Distance_To_Surface), '.eps');

        % Plot sliprate on X-T
        [Tsteps, Xsteps] = meshgrid(1e6 * time, 1e3 * x_up);
        h = pcolor(Xsteps', Tsteps', log10(SlipRateAtDist)');
        shading interp;

        if VitoColorFlag == 1
            plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofLogSlipRate_whole_Vito_surface', num2str(Distance_To_Surface), '.eps');
            colormap(flipud(black_rainbow_shear_long));
        end
        hold on;
        xline(VSregion(1), 'r' ,'linewidth', 2.0);
        xline(VSregion(2), 'r' ,'linewidth', 2.0);
        text(VSregion(1)+ 8, 38, 'Region 2', 'color', 'r', 'Fontsize', fontsize);
        
        % Add x_3
        text(-110, 100, strcat("$x_3$ = ", string(1000 * Distance_To_Surface - 5), "$\ \mathrm{mm}$"), 'color', 'w', 'Fontsize', fontsize);
        % Add the wave speeds
        
        % Several wave speed to show on the aplot
        cp = 2662.4;
        cs = 1279;
        nu = 0.35;
        cr = (0.874 + 0.196 * nu - 0.043 * nu^2 - 0.055 * nu^3) * cs;
        cX = [60, 80] - 15;
        cX = cX - norm(VSstart - WirePos1, 2) * 1e3;
        crY = [10, (cX(2) - cX(1)) * 1e3 / cr + 10];
        csY = [10, (cX(2) - cX(1)) * 1e3 / cs + 10];
        cpY = [10, (cX(2) - cX(1)) * 1e3 / cp + 10];
        
        % Add the wave speeds
        plot(cX, crY, 'w', 'linewidth', 2.0);
        text(cX(2) + 1, crY(2)+5, strcat('$c_r$ = 1.20 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
        plot(cX, csY, 'w', 'linewidth', 2.0);
        text(cX(2) + 1, crY(2) - 1, strcat('$c_s$ = 1.28 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
        plot(cX, cpY, 'w', 'linewidth', 2.0);
        text(cX(2) + 1, cpY(2) - 2, strcat('$c_p$ = 2.66 [km/s]'), 'color', 'w', 'Fontsize', fontsize - 5, 'interpreter', 'latex');
        hold off;
        
        set(h, 'EdgeColor', 'None');
        c = colorbar;
        clim([-15, 2]);
        ylabel(c,'Log of Slip rate [m/s]','FontName','Avenir','FontSize',fontsize, 'interpreter', 'latex');
        xlim(Xrange);
        ylim(Trange);
        title('X-T diagram of log slip rate')
        xlabel('$x_1\ \mathrm{[mm]}$', 'interpreter', 'latex');
        ylabel('Time [$\mathrm{\mu s}$]', 'interpreter', 'latex');
        % title(['X-T diagram of Slip rate at Depth = ', ' ', num2str(5 - 1e3 * Distance_To_Surface, '%.0f'), ' ', '[mm]'], 'interpreter', 'latex');
        set(gca, 'FontSize', fontsize);
        % daspect([yToxRatio, 1, 1]);
        disp(strcat('V_initial at $x_1 = 26$ mm: ', string(log10(SlipRateAtDist(1200, 2)))))
        
        % Save the figure
        if fmt == 'png'
            plotname = strcat(pwd, '/../PaperPlots_interp_1209/', videoprefix, '_X-TofLogSlipRate_whole_Vito_surface', num2str(Distance_To_Surface), '.png');
            print(figure(figNo), plotname, '-dpng', '-r500');
        else
            print('-vector', figure(figNo) ,plotname, '-depsc', '-r500');
        end
    end

  disp(strcat('V_initial at $x_1 = 26$ mm: ', string(log10(SlipRateAtDist(1200, 2)))));
end