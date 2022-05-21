function fV_vs_slip_time_function(totalprefix, target_x, z_location)
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
    pos(4) = pos(4) * 2;
    fig2.Position = pos;

    for i = 1:1:3
        videoprefix = strcat(num2str(i), totalprefix); 
        faultFileName = strcat('../faultFiles/', videoprefix, '-fault.h5');
        frontsurfFile = strcat('../frontsurfFiles/', videoprefix, '-frontsurf.h5');

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
        connection = h5read(faultFileName, '/topology/cells');
        connection = connection + 1;

        % Read nodal slip
        % Input the first wire position
        WirePos1 = [-0.025657; -0.014222; 0];

        % Calculate NodalXYZ2D 
        NodalXYZ2D = zeros(2, size(nodalXYZ, 2));
        for ii = 1:1:size(nodalXYZ, 2)
            NodalXYZ2D(1, ii) = sign(nodalXYZ(1, ii) - WirePos1(1)) * norm(nodalXYZ(1:2, ii) - WirePos1(1:2), 2);
            NodalXYZ2D(2, ii) = nodalXYZ(3, ii);
        end

        % Fault Range
        xrange = [-100, 150];
        
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
            end
        end
        
        
        %% Plot V-Slip history at X [mm]
        figure(1);
        for ii = 1:1:size(target_x, 2)
            % Generate sliprate-slip plot
            % Plot x and y axis ranges
            xrange = [0, 150];
            yrange = [0, 6];
            subplot(2, size(target_x, 2), size(target_x, 2) + ii);
            plot(QuerySlip(ii, :) * 1e6, QueryV(ii, :), 'linewidth', 2.0);
            hold on; grid on;
            scatter(QuerySlip(ii, :) * 1e6, QueryV(ii, :), 'filled');
            xlabel('Slip [$\mu$m]', 'interpreter', 'latex');
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
            scatter(QuerySlip(ii, :) * 1e6, -QueryShearStress(ii, :) ./ QueryNormalStress(ii, :), 'filled');

            if (ii == 1)
                ylabel('Friction coefficient', 'interpreter', 'latex');
            end
            xlim(xrange);
            ylim([0, 1]);
            probeLabel = strcat('X ={ }', num2str(target_x(ii)), '{ }[mm]');
            title(probeLabel, 'Fontsize', fontsize, 'interpreter', 'latex');
            set(gca, 'FontSize', fontsize);
        end


        %% -------------------------------------------------------------------------
        % Plot friction-sliprate time history at X [mm]
        % target_x same as above
        % Generate sliprate-time plot

        figure(2);
        xrange = [0, 110];
        % Generate sliprate-slip plot
        for ii = 1:1:size(target_x, 2)
            % Plot friction coefficient vs time
            subplot(2, size(target_x, 2), ii);
            plot(time * 1e6 - 10, -QueryShearStress(ii, :) ./ QueryNormalStress(ii, :), 'linewidth', 2.0);
            hold on; grid on;
            scatter(time * 1e6 - 10, -QueryShearStress(ii, :) ./ QueryNormalStress(ii, :), 'filled');

            if ii == 1
                ylabel('Friction coefficient', 'interpreter', 'latex');
            end

            xlim(xrange);
            ylim([0, 1]);
            probeLabel = strcat('X ={ }', num2str(target_x(ii)), '{ }[mm]');
            title(probeLabel, 'Fontsize', fontsize, 'interpreter', 'latex');
            set(gca, 'FontSize', fontsize);

            % Plot sliprate vs time
            subplot(2, size(target_x, 2), size(target_x, 2) + ii);
            plot(time * 1e6 - 10, QueryV(ii, :), 'linewidth', 2.0);
            hold on; grid on;
            scatter(time * 1e6 - 10, QueryV(ii, :), 'filled');
            xlabel('Time [$\mu$s]', 'interpreter', 'latex');
            if (ii == 1)
                ylabel('Slip rate [m/s]', 'interpreter', 'latex');
            end
            xlim(xrange);
            ylim(yrange);
            set(gca, 'FontSize', fontsize);
        end
    end
    figure(1);
    subplot(2, size(target_x, 2), 2 * size(target_x, 2));
    kids = get(gca, 'children');
    legend([kids(1), kids(3), kids(5)], 'mesh 3', 'mesh 2', 'mesh 1', 'location', 'best', 'interpreter', 'latex');
    set(gca, 'FontSize', fontsize);

    figure(2);
    subplot(2, size(target_x, 2), 2 * size(target_x, 2));
    kids = get(gca, 'children');
    legend([kids(1), kids(3), kids(5)], 'mesh 3', 'mesh 2', 'mesh 1', 'location', 'best', 'interpreter', 'latex');
    set(gca, 'FontSize', fontsize);

    % Save the files
    % plotname = strcat(pwd, '/../Vitoplots/', totalprefix, num2str(target_x), 'z_', num2str(z_location), '_fandV_VS_Slip.png');
    plotname = strcat(pwd, '/../Vitoplots/', totalprefix, 'z_', num2str(z_location), '_fandV_VS_Slip.png');
    disp(plotname);
    print(figure(1) ,plotname, '-dpng', '-r500');

    % plotname = strcat(pwd, '/../Vitoplots/', totalprefix, num2str(target_x),  'z_', num2str(z_location), '_fandV_VS_Time.png');
    plotname = strcat(pwd, '/../Vitoplots/', totalprefix, 'z_', num2str(z_location), '_fandV_VS_Time.png');
    disp(plotname);
    print(figure(2) ,plotname, '-dpng', '-r500');
end
