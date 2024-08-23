% Read results from hdf5 files.
setEnvironment;
% totalprefix = 'WithWallDRS1.5_1.5ModA0.008AmB0.005Load5_Vw2_fw0.1_theta0.036_-11_NULoad2dir0';
totalprefix = "W8_4ExcCorr2NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8"; 
% totalprefix = "W8_4ExcCorr2NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load8_Vw2_fw0.1_theta0.3399_-8_NULoad2dir0_duration200_8";
% totalprefix = "W8_4ExcCorr2NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load8_Vw2_fw0.1_theta0.8061_-9_NULoad2dir0_duration200_8";
% Node positions
target_x = [-100];

% Specify figure positions
fig1 = figure(1);
pos = fig1.Position;
pos(3) = pos(3) * size(target_x, 2);
pos(4) = pos(4) * 2;
fig1.Position = pos;

fig2 = figure(2);
pos = fig2.Position;
pos(3) = pos(3) * size(target_x, 2);
pos(4) = pos(4) * 4;
fig2.Position = pos;

for i = 2:1:2
    videoprefix = totalprefix; % strcat(num2str(i), totalprefix); 
    faultFileName = strcat('../faultFiles/', videoprefix, '-fault.h5');
    frontsurfFile = strcat('../frontsurfFiles/', videoprefix, '-frontsurf.h5');

    h5disp(faultFileName);
    h5disp(frontsurfFile);
    fontsize = 25;

    % Read time
    time = h5read(faultFileName, '/time');
    time = reshape(time, [1, size(time, 3)]);
    nOfTimeSteps = size(time, 2);

    % Read node geometry
    nodalXYZ = h5read(faultFileName, '/geometry/vertices');

    % Select the surface nodes
    I1 = (abs(nodalXYZ(3, :) - 0.005) <= 1e-7);
    nOfNodes = size(nodalXYZ, 2);

    % Read nodal slip, slip rate and traction
    Slip = h5read(faultFileName, '/vertex_fields/slip');
    SlipRate = h5read(faultFileName, '/vertex_fields/slip_rate');
    Traction = h5read(faultFileName, '/vertex_fields/traction');
    theta = h5read(faultFileName, '/vertex_fields/state_variable');
    % Connection = h5read(faultFileName, '/topology/cells');

    % Read nodal slip
    wireLabel = strcat('X = 0 [mm], wire position');
    VSLabel = 'VS Region [mm]';

    % Input the first wire position
    WirePos1 = [-0.025657; -0.014222; 0];

    VSstart = [0.006354, 0.003522, 0]';
    VSend = [0.063204, 0.035034, 0]';
    VSregion = 1e3 * [norm(VSstart - WirePos1, 2), norm(VSend - WirePos1, 2)];

    nOf2DNodes = sum(I1);

    FaultX = zeros(1, nOf2DNodes);
    surfaceNodesXYZ = nodalXYZ(:, I1);
    [~, I] = sort(surfaceNodesXYZ(1,:));

    % Plot x and y axis ranges
    xrange = [0, 150];
    yrange = [0, 20];
    yLogRange = [-10, 2];

    for i = 1:1:nOf2DNodes
        FaultX(i) = norm(surfaceNodesXYZ(1:2, i) - WirePos1(1:2), 2) * sign(surfaceNodesXYZ(1, i) - WirePos1(1));
    end
    FaultX = FaultX(I);
    surfaceNodesXYZ = surfaceNodesXYZ(1:2, I);

    % 2D node coordinates
    nodalXYZ2D = zeros(2, nOfNodes);
    for i = 1:1:nOfNodes
        nodalXYZ2D(1, i) = sign(nodalXYZ(1, i) - WirePos1(1)) ...
            * norm(nodalXYZ(1:2, i) - WirePos1(1:2), 2);
        nodalXYZ2D(2, i) = nodalXYZ(3, i);
    end

    % Magnitude of slip rate
    surfaceSlipRateMag = zeros(nOf2DNodes, nOfTimeSteps);
    surfaceSlipRate = zeros(3, nOf2DNodes, nOfTimeSteps);

    % Magnitude of slip
    surfaceSlipMag = zeros(nOf2DNodes, nOfTimeSteps);
    surfaceTheta = zeros(nOf2DNodes, nOfTimeSteps); 
    surfaceSlip = zeros(3, nOf2DNodes, nOfTimeSteps);

    % Shear and Normal stress
    surfaceShearStress = zeros(nOf2DNodes, nOfTimeSteps);
    surfaceNormalStress = zeros(nOf2DNodes, nOfTimeSteps);

    for t = 1:1:nOfTimeSteps
        surfaceSlipRate(:, :, t) = SlipRate(:, I1, t);
        surfaceSlip(:, :, t) = Slip(:, I1, t);
        surfaceShearStress(:, t) = Traction(1, I1, t);
        surfaceNormalStress(:, t) = Traction(3, I1, t);
        surfaceTheta(:, t) = theta(1, I1, t);

        for i = 1:1:nOf2DNodes
            surfaceSlipRateMag(i, t) = norm(surfaceSlipRate(:, i, t));
            surfaceSlipMag(i, t) = norm(surfaceSlip(:, i, t));
        end
        surfaceSlipRateMag(:, t) = surfaceSlipRateMag(I, t);
        surfaceSlipMag(:, t) = surfaceSlipMag(I, t);
        surfaceShearStress(:, t) = surfaceShearStress(I, t);
        surfaceNormalStress(:, t) = surfaceNormalStress(I, t);
        surfaceTheta(:, t) = surfaceTheta(I, t);
    end

    %% Plot V-Slip history at X [mm]
    figure(1);
    for ii = 1:1:size(target_x, 2)
        % Find the target index
        [~, Ind] = min(abs(FaultX - target_x(ii) / 1e3));

        % Generate sliprate-slip plot
        subplot(2, size(target_x, 2), size(target_x, 2) + ii);
        plot(surfaceSlipMag(Ind, 2:end) * 1e6, surfaceSlipRateMag(Ind, 2:end), 'linewidth', 2.0);
        hold on; grid on;
        scatter(surfaceSlipMag(Ind, 2:end) * 1e6, surfaceSlipRateMag(Ind, 2:end), 'filled');
        xlabel('Slip [$\mathrm{\mu m}$]');
        if (ii == 1) 
            ylabel('Slip rate [m/s]');
        end
        xlim(xrange);
        ylim(yrange);
        set(gca, 'FontSize', fontsize);

        % Generate friction coefficient-slip plot
        subplot(2, size(target_x, 2), ii);
        plot(surfaceSlipMag(Ind, 2:end) * 1e6, - surfaceShearStress(Ind, 2:end) ./ surfaceNormalStress(Ind, 2:end), 'linewidth', 2.0);
        hold on; grid on;
        scatter(surfaceSlipMag(Ind, 2:end) * 1e6, - surfaceShearStress(Ind, 2:end) ./ surfaceNormalStress(Ind, 2:end), 'filled');

        if (ii == 1)
            ylabel('Friction coefficient');
        end
        xlim(xrange);
        ylim([0, 1]);
        probeLabel = strcat('X ={ }', num2str(target_x(ii)), '{ }[mm]');
        title(probeLabel, 'Fontsize', fontsize);
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
        % Find the target index
        [~, Ind] = min(abs(FaultX - target_x(ii) / 1e3));

        % Plot friction coefficient vs time
        subplot(4, size(target_x, 2), ii);
        plot(time * 1e6 - 10, - surfaceShearStress(Ind, 1:end) ./ surfaceNormalStress(Ind, 2), 'linewidth', 2.0);
        hold on; grid on;
        scatter(time * 1e6 - 10, - surfaceShearStress(Ind, 1:end) ./ surfaceNormalStress(Ind, 2), 'filled');

        if ii == 1
            ylabel('Friction coefficient');
        end

        xlim(xrange);
        ylim([0, 1]);
        probeLabel = strcat('X ={ }', num2str(target_x(ii)), '{ }[mm]');
        title(probeLabel, 'Fontsize', fontsize);
        set(gca, 'FontSize', fontsize);

        % Plot sliprate vs time
        subplot(4, size(target_x, 2), size(target_x, 2) + ii);
        plot(time * 1e6 - 10, surfaceSlipRateMag(Ind, 1:end), 'linewidth', 2.0);
        hold on; grid on;
        scatter(time * 1e6 - 10, surfaceSlipRateMag(Ind, 1:end), 'filled');
        xlabel('Time [$\mathrm{\mu s}$]');
        if (ii == 1)
            ylabel('Slip rate [m/s]');
        end
        xlim(xrange);
        ylim(yrange);
        set(gca, 'FontSize', fontsize);

        % Plot logsliprate vs time
        subplot(4, size(target_x, 2), size(target_x, 2) * 2 + ii);
        plot(time * 1e6 - 10, log10(surfaceSlipRateMag(Ind, 1:end)), 'linewidth', 2.0);
        hold on; grid on;
        scatter(time * 1e6 - 10, log10(surfaceSlipRateMag(Ind, 1:end)), 'filled');
        xlabel('Time [$\mathrm{\mu s}$]');
        if (ii == 1)
            ylabel('Log slip rate [m/s]');
        end
        xlim(xrange);
        ylim(yLogRange);
        set(gca, 'FontSize', fontsize);
        
        % Plot slip vs time
        subplot(4, size(target_x, 2), size(target_x, 2) * 3 + ii);
        plot(time * 1e6 - 10, 1e6 * surfaceSlipMag(Ind, 1:end), 'linewidth', 2.0);
        hold on; grid on;
        scatter(time * 1e6 - 10, 1e6 * surfaceSlipMag(Ind, 1:end), 'filled');
        xlabel('Time [$\mathrm{\mu s}$]');
        if (ii == 1)
            ylabel('Slip [$\mathrm{\mu m}$]');
        end
        xlim(xrange);
        ylim([0, 100]);
        set(gca, 'FontSize', fontsize);

        % Find difference between f_peak and a * log(V/V_ini)
        fric = - surfaceShearStress(Ind, 1:end) ./ surfaceNormalStress(Ind, 2);
        V_this = surfaceSlipRateMag(Ind, 1:end);
        slip_this = surfaceSlipMag(Ind, 1:end);
        theta_this = surfaceTheta(Ind, 1:end);

        [val, ind] = max(fric);
        hypo_fric = fric(2) + 0.003 * log(V_this(ind) / V_this(2));
        disp("True, Hypothesized jump in friction coefficient: " + num2str(val - fric(2)) + ", " + num2str(hypo_fric - fric(2)));
        disp("slip, slip rate, theta initially: " + num2str(slip_this(2)) + ", " + num2str(V_this(2)) + ", " + num2str(theta_this(2)));
        disp("slip, slip rate, theta then: " + num2str(slip_this(ind)) + ", " + num2str(V_this(ind)) + ", " + num2str(theta_this(ind)));
    end
end
figure(1);
subplot(2, size(target_x, 2), 2 * size(target_x, 2));
kids = get(gca, 'children');
% legend([kids(1), kids(3), kids(5)], 'mesh 3', 'mesh 2', 'mesh 1', 'location', 'best');
set(gca, 'FontSize', fontsize);
set(gcf, 'color', 'w');

figure(2);
subplot(4, size(target_x, 2), 3 * size(target_x, 2));
kids = get(gca, 'children');
% legend([kids(1), kids(3), kids(5)], 'mesh 3', 'mesh 2', 'mesh 1', 'location', 'best');
set(gca, 'FontSize', fontsize);
set(gcf, 'color', 'w');

% Save the files
plotname = strcat(pwd, '/../Vitoplots/', totalprefix, '_fandV_VS_Slip.png');
print(figure(1) ,plotname, '-dpng', '-r500');

plotname = strcat(pwd, '/../Vitoplots/', totalprefix, '_fandV_VS_Time.png');
print(figure(2) ,plotname, '-dpng', '-r500');


