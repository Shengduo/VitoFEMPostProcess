% Calculate the xyz locations for inserting xyz values
setEnvironment;
load('../matFiles/THD58_Friction.mat', 'position_THD58_1');
load('../matFiles/realNormalStress.mat');
fontsize = 25;

% Set the plotting positions
% locs = 4:3:30; % For the gouge
locs = -30 : 5 : -10;  % For Homalite
time_slipPlt = 1e-6 * (0:10:190);  % For slip vs. time plot
% Set distance to surface
Distance_To_Surface = [0, 0.005];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% videoprefix = "2OPDirWithWallDRS1.5_1.5ModA0.011AmB0.005Load5_Vw2_fw0.1_theta0.006_-9_NULoad2dir0_duration200_0";
% videoprefix = "W16_4ExcCorr1OPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load10_Vw2_fw0.1_theta1.4664_-7_NULoad2dir0_duration200_10";
% videoprefix = "W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8";
videoprefix = "W8_4ExcCorr1NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2e+16_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8";

faultFileName = strcat('../faultFiles/', videoprefix, '-fault.h5');
% h5disp(faultFileName);
fontsize = 25;

% Read time
time = h5read(faultFileName, '/time');
time = reshape(time, [1, size(time, 3)]);
time = time - 10e-6;

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
videoprefix = strcat(videoprefix, "_Distance", num2str(Distance_To_Surface));

% Read nodal slip, slip rate and traction
SlipRate = h5read(faultFileName, '/vertex_fields/slip_rate');
Slip = h5read(faultFileName, '/vertex_fields/slip');
connection = h5read(faultFileName, '/topology/cells');
connection = connection + 1;
traction = h5read(faultFileName, '/vertex_fields/traction');

SlipRate = SlipRate(:, :, 2:end);
traction = traction(:, :, 2:end);
Slip = Slip(:, :, 2:end);
time = time(2:end);


nOfTimeSteps = size(time, 2);

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
slipMag = zeros(nOfNodes, nOfTimeSteps);
slipRateMag = zeros(nOfNodes, nOfTimeSteps);

for t = 1:1:nOfTimeSteps
    for i = 1:1:nOfNodes
        slipRateMag(i, t) = norm(SlipRate(:, i, t), 2);
        slipMag(i, t) = norm(Slip(:, i, t), 2);
    end
end


% Get the mesh stored as triangulation
TR = triangulation(connection', NodalXYZ2D');

%% Get the interpolation of slip at these locations
% Get the start and end location of  VS region
VS_start_x = norm(VSstart - VSstart, 2);
VS_end_x = norm(VSend - VSstart, 2);
Fault_start_x = -norm(FaultStart - VSstart, 2);
Fault_end_x = norm(FaultEnd - VSstart, 2);

% Camera specific values
Pix_10mm = 81;  %Number of pixels in a grid of 10 mm ;  144 pixels in 2 mm -> 72 pixels in 1 mm
pxsize=10/(Pix_10mm)*1e-3; % pixel size in m/pixel  (Pix_10mm pixels  in 10 mm)
stepsize = 1;

% First specific the x and y's
x_up = locs * 1e-3;

% x_up for slip plots
x_up_slipPlt = min(NodalXYZ2D(1, :)) : 1e-3 : max(NodalXYZ2D(1, :));
x_slipPlt_range = [min(NodalXYZ2D(1, :)), max(NodalXYZ2D(1, :))];
slipPlt_range = [0, 350];

% Generate the first plot
figNo = 1;
fig1 = figure(figNo);
fig1.Position = [73 1 2488 1321];
figNo = figNo + 1;

fig2 = figure(figNo);
fig2.Position = [73, 1, 2488, 1321];
for j = 1:1:length(Distance_To_Surface)
    % Get the querying positions
    [Xq, Yq] = meshgrid(x_up, Distance_To_Surface(j));
    P = [Xq', Yq'];
    
    [Xq_slipPlt, Yq_slipPlt] = meshgrid(x_up_slipPlt, Distance_To_Surface(j));
    P_slipPlt = [Xq_slipPlt', Yq_slipPlt'];

    elementID = pointLocation(TR, P);
    elementID_slipPlt = pointLocation(TR, P_slipPlt);

    SlipRateAtDist = zeros(size(Xq, 2), size(time, 2));
    SlipAtDist = zeros(size(Xq, 2), size(time, 2));
    SlipAtDist_slipPlt = zeros(size(Xq_slipPlt, 2), size(time, 2));
    ShearStressAtDist = zeros(size(Xq, 2), size(time, 2));
    for t = 1:1:size(time, 2)
        for ele = 1:1:size(elementID, 1)
            vel_x = scatteredInterpolant(NodalXYZ2D(:, connection(:, elementID(ele)))', squeeze(slipRateMag(connection(:, elementID(ele)), t)), 'natural');
            slip_x = scatteredInterpolant(NodalXYZ2D(:, connection(:, elementID(ele)))', squeeze(slipMag(connection(:, elementID(ele)), t)), 'natural');
            stress_x = scatteredInterpolant(NodalXYZ2D(:, connection(:, elementID(ele)))', squeeze(traction(1, connection(:, elementID(ele)), t))', 'natural');
            SlipRateAtDist(ele, t) = vel_x(x_up(ele), Distance_To_Surface(j));
            SlipAtDist(ele, t) = slip_x(x_up(ele), Distance_To_Surface(j));
            ShearStressAtDist(ele, t) = stress_x(x_up(ele), Distance_To_Surface(j));
        end
    end

    for t = 1:1:size(time, 2)
        for ele = 1:1:size(elementID_slipPlt, 1)
            slip_x_slipPlt = scatteredInterpolant(NodalXYZ2D(:, connection(:, elementID_slipPlt(ele)))', squeeze(slipMag(connection(:, elementID_slipPlt(ele)), t)), 'natural');
            SlipAtDist_slipPlt(ele, t) = slip_x_slipPlt(x_up_slipPlt(ele), Distance_To_Surface(j));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Define ranges for the plots
    slipRange = [0, 40];
    slipRateRange = [0, 30];
    frictionRange = [0, 1];
    
    % Compute si0, tau0 [MPa]
    P = 14.3;
    si0 = P * cosd(29)^2;
    tau0 = P * cosd(29) * sind(29);
    
    % Compute change of shear stress
    shearChange_seq1 = 1e-6 * (ShearStressAtDist - ShearStressAtDist(:, 1));
    shearChangeRange = 1.2 * [min(min(shearChange_seq1)), max(max(shearChange_seq1))];
    
    % Compute dtau / Ninferred
    si_interp = interp1(position_THD58_1, si(1, :), locs);
    dTauBySi = shearChange_seq1 ./ si_interp(1, :)';
    % dTauBySi = shearChange_seq1 ./ si0;
    
    dTauBySiRange = 1.2 * [min(min(dTauBySi)), max(max(dTauBySi))];

    % Plot and save tauChange vs. sliprate and slip
    for i = 1:1:length(locs)
        figure(1);
        %% delta tau vs. sliprate and slip, seq 1
        % Adjust the figure size
        subplot(2, length(locs), i + (j - 1) * length(locs))
    
        % \delta tau vs. slip
        plot(1e6 * SlipAtDist(i, :), shearChange_seq1(i, :), 'linewidth', 2.0);
        hold on; grid on;
        scatter(1e6 * SlipAtDist(i, :), shearChange_seq1(i, :), 'filled');
        
        if j == length(Distance_To_Surface)
            xlabel('Slip [$\mathrm{\mu m}$]');
        end

        if j == 1
            title(strcat('$x_1$ = ', num2str(locs(i))), 'interpreter', 'latex');
        end

        if i == 1
            ylabel('$\Delta \tau$ [MPa]'); 
        end

        xlim(slipRange);
        ylim([-4, 4]);
        set(gca, 'fontsize', fontsize);
    
%         % \delta \tau vs. slip rate
%         subplot(2, 1, 2);
%         plot(SlipRateAtDist(i, :), shearChange_seq1(i, :), 'linewidth', 2.0);
%         hold on; grid on;
%         scatter(SlipRateAtDist(i, :), shearChange_seq1(i, :), 'filled');
%         xlabel('Slip rate [m/s]');
%         ylabel('$\Delta \tau$ [MPa]');
%         % title(['Slip rate vs. slip at x = ', num2str(position_THD58_1(idx)), ' [mm] (Seq 1)']);
%         xlim(slipRateRange);
%         ylim(shearChangeRange);
%         set(gca, 'fontsize', fontsize);
%         savename = "../plots20221117/UniformSlipWeakening/DTauAt" + num2str(locs(i)) + ".png";
%         print(fig, savename, '-dpng', '-r500');
%         figNo = figNo + 1;
        
%         %% d tau / si inferred vs slip rate, seq 1
%         % Adjust the figure size
%         fig = figure(figNo);
%         pos = fig.Position;
%         pos(4) = pos(4) * 3.0;
%         pos(3) = pos(3) * 1.5;
%         fig.Position = pos;
%         
%         
%         % vs slip
%         subplot(2, 1, 1);
%         plot(1e6 * SlipAtDist(i, :), dTauBySi(i, :), 'linewidth', 2.0);
%         hold on; grid on;
%         scatter(1e6 * SlipAtDist(i, :), dTauBySi(i, :), 'filled');
%         xlabel('Slip [$\mathrm{\mu m}$]');
%         % ylabel('$\Delta \tau / \sigma_{inferred}$');
%         % title(['Real $\Delta \tau / \sigma_{inferred}$  at x = ', num2str(locs(i)), ' [mm]']);
%         ylabel('$\Delta \tau / \sigma_{0}$');
%         title(['Real $\Delta \tau / \sigma_{0}$  at x = ', num2str(locs(i)), ' [mm]']);
%         xlim(slipRange);
%         % ylim(dTauBySiRange);
%         ylim([-4, 4]);
%         set(gca, 'fontsize', fontsize);
        
%         % vs slip rate
%         subplot(2, 1, 2);
%         plot(SlipRateAtDist(i, :), dTauBySi(i, :), 'linewidth', 2.0);
%         hold on; grid on;
%         scatter(SlipRateAtDist(i, :), dTauBySi(i, :), 'filled');
%         xlabel('Slip rate [m/s]');
%         ylabel('$\Delta \tau / \sigma_{inferred}$');
%         % ylabel('$\Delta \tau / \sigma_{0}$');
%         % title(['Exp $$ vs. slip rate at x = ', num2str(target_x), ' [mm] (Seq 1)']);
%         xlim(slipRateRange);
%         ylim(dTauBySiRange);
%         set(gca, 'fontsize', fontsize);
%         savename = "../plots20221117/UniformSlipWeakening/DTauByNAt" + num2str(locs(i)) + ".png";
%         print(fig, savename, '-dpng', '-r500');
    
%         figNo = figNo + 1;
%         close all;
    end

% #############################################################################
% Displacement vs.location plot
    figure(2);
%    lines = [];
    colors = flipud(autumn(length(time_slipPlt)));
    slipAtTimes = interp1(time, SlipAtDist_slipPlt', time_slipPlt);
    for k = 1:1:length(time_slipPlt)
        subplot(2, 1, j);
        ln = plot(1e3 * x_up_slipPlt, 1e6 * slipAtTimes(k, :), 'LineWidth', 2.0, 'Color', colors(k, :));
%         lines = [lines, ln];
        hold on;
        grid on;
        ylabel('Slip $[\mathrm{\mu m}]$', 'FontSize', 25);
        if j == length(Distance_To_Surface)
            xlabel('$x_1\ [\mathrm{mm}]$', 'FontSize', 25);
        end
        xlim(1e3 * x_slipPlt_range);
        ylim(slipPlt_range);
        set(gca, 'fontsize', fontsize);
    end

    % add legend
    Lgnd = legend(strcat('$\ t =\ $', string(num2cell(1e6 * time_slipPlt)), '$\ \mathrm{\mu s}$'), 'interpreter', 'latex');
    set(Lgnd,...
    'Position',[0.909179534830337 0.253744488954795 0.0891623665665501 0.499242997728993]);
end

% Give entire title
figure(1);
sgtitle("$V_{ini}=10^{-7}\ \mathrm{m/s}$, $x_3 = -5\ [\mathrm{mm}]$ (up row) and $x_3 = 0\ [\mathrm{mm}]$ (lower row)", 'fontsize', 25);
set(gcf, 'color', 'w');
plotPath = "/home/shengduo/Desktop/photos/" + videoprefix + "-1.png";
print(figure(1), plotPath, '-dpng', '-r300');

% Give entire title
figure(2);
sgtitle("$V_{ini}=10^{-7}\ \mathrm{m/s}$, $x_3 = -5\ [\mathrm{mm}]$ (up row) and $x_3 = 0\ [\mathrm{mm}]$ (lower row)", 'fontsize', 25);
set(gcf, 'color', 'w');
plotPath = "/home/shengduo/Desktop/photos/" + videoprefix + "-2.png";
print(figure(2), plotPath, '-dpng', '-r300');

close all;
