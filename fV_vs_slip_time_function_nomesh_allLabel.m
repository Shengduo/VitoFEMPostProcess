function fV_vs_slip_time_function_nomesh_allLabel(totalprefixs, target_x, z_location, theta_flag, legendFlag)
    % Read results from hdf5 files.
    % totalprefix = 'WithWallDRS1.5_1.5ModA0.008AmB0.005Load5_Vw2_fw0.1_theta0.036_-11_NULoad2dir0';

    % Node positions
    % target_x = [-50, -10, -5];

    % Specify figure positions
    
    fig1 = figure(1);
    pos = fig1.Position;
    pos(3) = pos(3) * size(target_x, 2);
    pos(4) = pos(4) * 2;
    fig1.Position = pos;

    fig2 = figure(2);
    pos = fig2.Position;
    pos(3) = pos(3) * size(target_x, 2);
    if theta_flag == 1
        pos(4) = pos(4) * 4;
    else
        pos(4) = pos(4) * 3;
    end
    fig2.Position = pos;
    
    % Store fracture energy
    high_f = zeros(1, length(target_x)); 
    mean_high_f = 0.;
    
    low_f = zeros(1, length(target_x)); 
    mean_low_f = 0.;

    f_energy = zeros(1, length(target_x));
    mean_f_energy = 0.;

    for i = 1:1:length(totalprefixs)
        % videoprefix = strcat(num2str(i), totalprefix); 
        videoprefix = totalprefixs(i);
        faultFileName = strcat('../faultFiles/', videoprefix, '-fault.h5');

        % h5disp(faultFileName);
        % h5disp(frontsurfFile);
        fontsize = 25;

        % Read time
        time = h5read(faultFileName, '/time');
        time = reshape(time, [1, size(time, 3)]);
        nOfTimeSteps = size(time, 2);

        % Read node geometry
        nodalXYZ = h5read(faultFileName, '/geometry/vertices');
        nOfNodes = size(nodalXYZ, 2);
        
        % Read nodal slip, slip rate and traction
        slip = h5read(faultFileName, '/vertex_fields/slip');
        SlipRate = h5read(faultFileName, '/vertex_fields/slip_rate');
        traction = h5read(faultFileName, '/vertex_fields/traction');
        if theta_flag == 1
            theta = h5read(faultFileName, '/vertex_fields/state_variable');
        end
        
        connection = h5read(faultFileName, '/topology/cells');
        connection = connection + 1;

        % Read nodal slip
        % Input the first wire position
        WirePos1 = [-0.025657; -0.014222; 0];
        VSstart = [0.006354, 0.003522, 0]';
        % Calculate NodalXYZ2D 
        NodalXYZ2D = zeros(2, size(nodalXYZ, 2));
        for ii = 1:1:size(nodalXYZ, 2)
            NodalXYZ2D(1, ii) = sign(nodalXYZ(1, ii) - VSstart(1)) * norm(nodalXYZ(1:2, ii) - VSstart(1:2), 2);
            NodalXYZ2D(2, ii) = nodalXYZ(3, ii);
        end

        % Fault Range
        xrange = [-100, 150] - 1e3 * norm(VSstart - WirePos1);
        
        % Magnitude of slip
        slipMag = zeros(nOfNodes, nOfTimeSteps);
        
        % Magnitude of slip rate
        slipRateMag = zeros(nOfNodes, nOfTimeSteps);

        for t = 1:1:nOfTimeSteps
            for ii = 1:1:nOfNodes
                slipRateMag(ii, t) = norm(SlipRate(:, ii, t), 2);
                slipMag(ii, t) = norm(slip(:, ii, t), 2);
            end
        end

        % Get the mesh stored as triangulation, get the location of places
        % of interest
        TR = triangulation(connection', NodalXYZ2D');
        QueryXYZ = [target_x / 1e3; z_location / 1e3 * ones(1, size(target_x, 2))]';
        elementID = pointLocation(TR, QueryXYZ);
        nOfQueryPts = size(target_x, 2);
        
        %% Get sliprate, shear stress aned normal stress history at target points
        QuerySlip = zeros(nOfQueryPts, nOfTimeSteps);
        QueryV = zeros(nOfQueryPts, nOfTimeSteps);
        QueryShearStress = ones(nOfQueryPts, nOfTimeSteps);
        QueryNormalStress = ones(nOfQueryPts, nOfTimeSteps);
        if theta_flag == 1
            QueryTheta = ones(nOfQueryPts, nOfTimeSteps);
        end
        
        % Get the interpolation
        for t = 1:1:nOfTimeSteps
            for pt = 1:1:nOfQueryPts
                % Slip
                slip = scatteredInterpolant(NodalXYZ2D(:, connection(:, elementID(pt)))', squeeze(slipMag(connection(:, elementID(pt)), t)), 'natural');
                QuerySlip(pt, t) = slip(QueryXYZ(pt, 1), QueryXYZ(pt, 2));
                
                % Slip rate
                vel_x = scatteredInterpolant(NodalXYZ2D(:, connection(:, elementID(pt)))', squeeze(slipRateMag(connection(:, elementID(pt)), t)), 'natural');
                QueryV(pt, t) = vel_x(QueryXYZ(pt, 1), QueryXYZ(pt, 2));
                
                % Shear stress
                shear = scatteredInterpolant(NodalXYZ2D(:, connection(:, elementID(pt)))', squeeze(traction(1, connection(:, elementID(pt)), t))', 'natural');
                QueryShearStress(pt, t) = shear(QueryXYZ(pt, 1), QueryXYZ(pt, 2));
                
                % Normal stress
                normal = scatteredInterpolant(NodalXYZ2D(:, connection(:, elementID(pt)))', squeeze(traction(3, connection(:, elementID(pt)), t))', 'natural');
                QueryNormalStress(pt, t) = normal(QueryXYZ(pt, 1), QueryXYZ(pt, 2));
                
                if theta_flag == 1
                    % State variable
                    stateVar = scatteredInterpolant(NodalXYZ2D(:, connection(:, elementID(pt)))', squeeze(theta(1, connection(:, elementID(pt)), t))', 'natural');
                    QueryTheta(pt, t) = stateVar(QueryXYZ(pt, 1), QueryXYZ(pt, 2));
                end
            end
        end
        
        % Save QueryNormalStress, QueryV
        save('~/InverseProblems/RateStateFrictionTraining/data/realData.mat', 'QueryV', 'QueryNormalStress', 'time');
        
        %% Plot V-Slip history at X [mm]
        figure(1);
        for ii = 1:1:size(target_x, 2)
            % Generate sliprate-slip plot
            % Plot x and y axis ranges
            xrange = [0, 100];
            yrange = [0, 10];
            subplot(2, size(target_x, 2), size(target_x, 2) + ii);
            plot(QuerySlip(ii, :) * 1e6, QueryV(ii, :), 'linewidth', 2.0);
            hold on; grid on;
            if legendFlag ~= 2
                scatter(QuerySlip(ii, :) * 1e6, QueryV(ii, :), 'filled');
            end
            xlabel('Slip [$\mathrm{\mu m}$]', 'interpreter', 'latex');
            if (ii == 1) 
                ylabel('Slip rate [m/s]', 'interpreter', 'latex');
            end
            xlim(xrange);
            ylim(yrange);
            set(gca, 'FontSize', fontsize);

            % Generate friction coefficient-slip plot
            subplot(2, size(target_x, 2), ii);
            plot(QuerySlip(ii, :) * 1e6, -QueryShearStress(ii, :) ./ QueryNormalStress(ii, :), 'linewidth', 2.0);
            hold on; grid on;
            if legendFlag ~= 2
                scatter(QuerySlip(ii, :) * 1e6, -QueryShearStress(ii, :) ./ QueryNormalStress(ii, :), 'filled');
            end
            xlabel('Slip [$\mathrm{\mu m}$]', 'interpreter', 'latex');
            if (ii == 1)
                ylabel('Friction coefficient', 'interpreter', 'latex');
            end
            xlim(xrange);
            ylim([0, 1]);
            probeLabel = strcat('$x_1$ ={ }', num2str(target_x(ii)), '{ }[mm]');
            title(probeLabel, 'Fontsize', fontsize, 'interpreter', 'latex');
            set(gca, 'FontSize', fontsize);
            set(gcf,'color','white');

            % Calculate friction energy
            high_f(ii) = max(-QueryShearStress(ii, :) ./ QueryNormalStress(ii, :));
            low_f(ii) = min(-QueryShearStress(ii, :) ./ QueryNormalStress(ii, :));
            f_energy(ii) = trapz(QuerySlip(ii, 2:end) * 1e6, -QueryShearStress(ii, 2:end) ./ QueryNormalStress(ii, 2:end) - low_f(ii));
        end


        %% -------------------------------------------------------------------------
        % Plot friction-sliprate time history at X [mm]
        % target_x same as above
        % Generate sliprate-time plot

        figure(2);
        xrange = [0, 110];
        if theta_flag == 1
            nrows = 4;
        else
            nrows = 3;
        end
        % Generate sliprate-slip plot
        for ii = 1:1:size(target_x, 2)
            % Plot friction coefficient vs time
            subplot(nrows, size(target_x, 2), ii);
            plot(time * 1e6 - 10, -QueryShearStress(ii, :) ./ QueryNormalStress(ii, :), 'linewidth', 2.0);
            hold on; grid on;
            xlabel('Time [$\mathrm{\mu s}$]', 'interpreter', 'latex');
            if legendFlag ~= 2
                scatter(time * 1e6 - 10, -QueryShearStress(ii, :) ./ QueryNormalStress(ii, :), 'filled');
            end

            if ii == 1
                ylabel('Friction coefficient', 'interpreter', 'latex');
            end

            xlim(xrange);
            ylim([0, 1]);
            probeLabel = strcat('$x_1$ ={ }', num2str(target_x(ii)), '{ }[mm]');
            title(probeLabel, 'Fontsize', fontsize, 'interpreter', 'latex');
            set(gca, 'FontSize', fontsize);

            % Plot sliprate vs time
            subplot(nrows, size(target_x, 2), size(target_x, 2) + ii);
            plot(time * 1e6 - 10, QueryV(ii, :), 'linewidth', 2.0);
            hold on; grid on;
            if legendFlag ~= 2
                scatter(time * 1e6 - 10, QueryV(ii, :), 'filled');
            end
            xlabel('Time [$\mathrm{\mu s}$]', 'interpreter', 'latex');
            if (ii == 1)
                ylabel('Slip rate [m/s]', 'interpreter', 'latex');
            end
            xlim(xrange);
            ylim(yrange);
            set(gca, 'FontSize', fontsize);
            
            % Plot slip vs time
            subplot(nrows, size(target_x, 2), 2 * size(target_x, 2) + ii);
            plot(time * 1e6 - 10, 1e6 * QuerySlip(ii, :), 'linewidth', 2.0);
            hold on; grid on;
            if legendFlag ~= 2
                scatter(time * 1e6 - 10, 1e6 * QuerySlip(ii, :), 'filled');
            end
            if nrows == 3
                xlabel('Time [$\mathrm{\mu s}$]', 'interpreter', 'latex');
            end
            if (ii == 1)
                ylabel('Slip [$\mathrm{\mu m}$]', 'interpreter', 'latex');
            end
            xlim(xrange);
            ylim([0, 50]);
            set(gca, 'FontSize', fontsize);
            
            % Plot theta vs time
            if theta_flag == 1
                subplot(nrows, size(target_x, 2), 3 * size(target_x, 2) + ii);
                plot(time * 1e6 - 10, log10(QueryTheta(ii, :)), 'linewidth', 2.0);
                hold on; grid on;
                if legendFlag ~= 2
                    scatter(time * 1e6 - 10, log10(QueryTheta(ii, :)), 'filled');
                end
                xlabel('Time [$\mathrm{\mu s}$]', 'interpreter', 'latex');
                if (ii == 1)
                    ylabel('$\mathrm{\log(\theta)}$ [$\mathrm{s}$]', 'interpreter', 'latex');
                end
                xlim(xrange);
                ylim([-10, 10]);
                set(gca, 'FontSize', fontsize);
            end
            
        end
        
        
    end
    
    %% Shit 
    load('../matFiles/2DFricQuery.mat');
    figure(1);
    subplot(2, 2, 1);
    plot(1e6 * QuerySlip(1, :), Queryfric(1, :), 'linewidth', 1.5);
    subplot(2, 2, 2);
    plot(1e6 * QuerySlip(2, :), Queryfric(2, :), 'linewidth', 1.5);
    legend('3D', '2D', 'fontsize', fontsize, 'location', 'best', 'interpreter', 'latex');

    figure(2);
    subplot(3, 2, 3);
    plot(QueryTime, QueryV(1, :), 'linewidth', 1.5);
    subplot(3, 2, 4);
    plot(QueryTime, QueryV(2, :), 'linewidth', 1.5);
    legend('3D', '2D', 'fontsize', fontsize, 'location', 'best', 'interpreter', 'latex');

    figure(1);
    if legendFlag == 1
        subplot(2, size(target_x, 2), 2 * size(target_x, 2));
        kids = get(gca, 'children');
        legend([kids(1), kids(3), kids(5)], 'mesh 3', 'mesh 2', 'mesh 1', 'location', 'best', 'interpreter', 'latex');
    end
    if legendFlag == 2
        subplot(2, size(target_x, 2), size(target_x, 2));
        kids = get(gca, 'children');
        % legend('Ageing', 'Slip', 'location', 'best', 'interpreter', 'latex');
    end
    set(gcf,'color','white');
    set(gca, 'FontSize', fontsize);

    figure(2);
    if legendFlag == 1
        subplot(3, size(target_x, 2), 3 * size(target_x, 2));
        kids = get(gca, 'children');
        legend([kids(1), kids(3), kids(5)], 'mesh 3', 'mesh 2', 'mesh 1', 'location', 'best', 'interpreter', 'latex');
    end
    if legendFlag == 2
        subplot(3, size(target_x, 2), 2 * size(target_x, 2));
        kids = get(gca, 'children');
        % legend('Ageing', 'Slip', 'location', 'best', 'interpreter', 'latex');
    end
    set(gcf,'color','white');
    set(gca, 'FontSize', fontsize);

    % Save the files
    % plotname = strcat(pwd, '/../Vitoplots/', totalprefix, num2str(target_x), 'z_', num2str(z_location), '_fandV_VS_Slip.eps');
    plotname = strcat(pwd, '/../fVSlipTimePlots/', totalprefixs(1), 'z_', num2str(z_location), '_fandV_VS_Slip.eps');
    disp(plotname);
    print('-vector', figure(1) ,plotname, '-depsc');

    % plotname = strcat(pwd, '/../Vitoplots/', totalprefix, num2str(target_x),  'z_', num2str(z_location), '_fandV_VS_Time.eps');
    plotname = strcat(pwd, '/../fVSlipTimePlots/', totalprefixs(1), 'z_', num2str(z_location), '_f_SlipandV_VS_Time.eps');
    if theta_flag == 1
        plotname = strcat(pwd, '/../fVSlipTimePlots/', totalprefixs(1), 'z_', num2str(z_location), '_f_Slip_VandTheta_VS_Timegpng');
    end
    disp(plotname);
    print('-vector', figure(2) ,plotname, '-depsc');

    %% Report fracture energy
    mean_low_f = mean(low_f);
    mean_high_f = mean(high_f);
    mean_f_energy = mean(f_energy);

    disp(strcat("high_f: ", num2str(high_f)));
    disp(strcat("mean_high_f: ", num2str(mean_high_f)));
    disp(strcat("low_f: ", num2str(low_f)));
    disp(strcat("mean_low_f: ", num2str(mean_low_f)));
    disp(strcat("f_energy: ", num2str(f_energy)));
    disp(strcat("mean_f_energy: ", num2str(mean_f_energy)));

    % Compute effective slipping distance D0 for linear slip weakening
    D0 = 2 * mean_f_energy / (mean_high_f - mean_low_f);
    disp(strcat("Effective slip weakening D0: ", num2str(D0)));


    
end
