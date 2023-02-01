% Calculate the xyz locations for inserting xyz values
setEnvironment;
load('../matFiles/THD58_Friction.mat');
load('../matFiles/realNormalStress.mat');
fontsize = 25;

% Define ranges for the plots
slipRange = [0, 100];
slipRateRange = [0, 30];
frictionRange = [0, 1];

% Compute si0, tau0 [MPa]
P = 14.3;
si0 = P * cosd(29)^2;
tau0 = P * cosd(29) * sind(29);

% Compute change of shear stress
shearChange_seq1 = friction_AVG_THD58_1 * si0 - tau0;
shearChangeRange = 1.2 * [min(min(shearChange_seq1)), max(max(shearChange_seq1))];

% Compute dtau / Ninferred
dTauBySi = shearChange_seq1 ./ si(1, :)';
dTauBySiRange = 1.2 * [min(min(dTauBySi)), max(max(dTauBySi))];

% Set the plotting positions
locs = 4:3:30;

% Plot and save tauChange vs. sliprate and slip
for i = 1:1:length(locs)
    % Find the index of location
    [~, idx] = min(abs(position_THD58_1 - locs(i))); 
    
    figNo = 1;
    %% delta tau vs. sliprate and slip, seq 1
    % Adjust the figure size
    fig = figure(figNo);
    pos = fig.Position;
    pos(4) = pos(4) * 3.0;
    pos(3) = pos(3) * 1.5;
    fig.Position = pos;

    % \delta tau vs. slip
    subplot(2, 1, 1);
    plot(Slip_avg_THD58_1(idx, :), shearChange_seq1(idx, :), 'linewidth', 2.0);
    hold on; grid on;
    scatter(Slip_avg_THD58_1(idx, :), shearChange_seq1(idx, :), 'filled');
    xlabel('Slip [$\mathrm{\mu m}$]');
    ylabel('$\Delta \tau$ [MPa]');
    title(['Exp $\Delta \tau$ at ', num2str(position_THD58_1(idx)), ' [mm]']);
    xlim(slipRange);
    ylim(shearChangeRange);
    set(gca, 'fontsize', fontsize);

    % \delta \tau vs. slip rate
    subplot(2, 1, 2);
    plot(Slip_rate_avg_THD58_1(idx, :), shearChange_seq1(idx, 2:end), 'linewidth', 2.0);
    hold on; grid on;
    scatter(Slip_rate_avg_THD58_1(idx, :), shearChange_seq1(idx, 2:end), 'filled');
    xlabel('Slip rate [m/s]');
    ylabel('$\Delta \tau$ [MPa]');
    % title(['Slip rate vs. slip at x = ', num2str(position_THD58_1(idx)), ' [mm] (Seq 1)']);
    xlim(slipRateRange);
    ylim(shearChangeRange);
    set(gca, 'fontsize', fontsize);
    savename = "../plots20221117/experiment/DTauAt" + num2str(locs(i)) + ".png";
    print(fig, savename, '-dpng', '-r500');
    figNo = figNo + 1;
    
    %% d tau / si inferred vs slip rate, seq 1
    % Adjust the figure size
    fig = figure(figNo);
    pos = fig.Position;
    pos(4) = pos(4) * 3.0;
    pos(3) = pos(3) * 1.5;
    fig.Position = pos;
    
    
    % vs slip
    subplot(2, 1, 1);
    plot(Slip_avg_THD58_1(idx, :), dTauBySi(idx, :), 'linewidth', 2.0);
    hold on; grid on;
    scatter(Slip_avg_THD58_1(idx, :), dTauBySi(idx, :), 'filled');
    xlabel('Slip [$\mathrm{\mu m}$]');
    ylabel('$\Delta \tau / \sigma_{inferred}$');
    title(['Exp $\Delta \tau / \sigma_{inferred}$  at x = ', num2str(position_THD58_1(idx)), ' [mm]']);
    xlim(slipRange);
    ylim(dTauBySiRange);
    set(gca, 'fontsize', fontsize);
    
    % vs slip rate
    subplot(2, 1, 2);
    plot(Slip_rate_avg_THD58_1(idx, :), dTauBySi(idx, 2:end), 'linewidth', 2.0);
    hold on; grid on;
    scatter(Slip_rate_avg_THD58_1(idx, :), dTauBySi(idx, 2:end), 'filled');
    xlabel('Slip rate [m/s]');
    ylabel('$\Delta \tau / \sigma_{inferred}$');
    % title(['Exp $$ vs. slip rate at x = ', num2str(target_x), ' [mm] (Seq 1)']);
    xlim(slipRateRange);
    ylim(dTauBySiRange);
    set(gca, 'fontsize', fontsize);
    savename = "../plots20221117/experiment/DTauByNAt" + num2str(locs(i)) + ".png";
    print(fig, savename, '-dpng', '-r500');


    figNo = figNo + 1;

    close all;
end

% %% f and sliprate vs. slip, seq 2
% % Adjust the figure size
% fig = figure(figNo);
% pos = fig.Position;
% pos(4) = pos(4) * 3;
% pos(3) = pos(3) * 1.5;
% fig.Position = pos;
% 
% % f vs. slip
% subplot(2, 1, 1);
% plot(Slip_avg_THD58_2(id, :), friction_AVG_THD58_2(id, :), 'linewidth', 2.0);
% hold on; grid on;
% scatter(Slip_avg_THD58_2(id, :), friction_AVG_THD58_2(id, :), 'filled');
% xlabel('Slip [$\mu$ m]');
% ylabel('Friction coefficient');
% title(['Friction vs. slip at x = ', num2str(target_x), ' [mm] (Seq 2)']);
% xlim(slipRange);
% ylim(frictionRange);
% set(gca, 'fontsize', fontsize);
% 
% % V vs. slip
% subplot(2, 1, 2);
% plot(Slip_avg_THD58_2(id, 2:end), Slip_rate_avg_THD58_2(id, :), 'linewidth', 2.0);
% hold on; grid on;
% scatter(Slip_avg_THD58_2(id, 2:end), Slip_rate_avg_THD58_2(id, :), 'filled');
% xlabel('Slip [$\mu$ m]');
% ylabel('Slip rate [m/s]');
% title(['Slip rate vs. slip at x = ', num2str(target_x), ' [mm] (Seq 2)']);
% xlim(slipRange);
% ylim(slipRateRange);
% set(gca, 'fontsize', fontsize);
% savename = "../Vitoplots/fV_slip_x" + num2str(target_x) + "seq2.png";
% print(fig, savename, '-dpng', '-r500');
% 
% 
% figNo = figNo + 1;
% %% f vs slip rate, seq 2
% % Adjust the figure size
% fig = figure(figNo);
% pos = fig.Position;
% pos(4) = pos(4) * 1.5;
% pos(3) = pos(3) * 1.5;
% fig.Position = pos;
% plot(Slip_rate_avg_THD58_2(id, :), friction_AVG_THD58_2(id, 2:end), 'linewidth', 2.0);
% hold on; grid on;
% scatter(Slip_rate_avg_THD58_2(id, :), friction_AVG_THD58_2(id, 2:end), 'filled');
% xlabel('Slip rate [m/s]');
% ylabel('Friction coefficient');
% title(['Friction vs. slip rate at x = ', num2str(target_x), ' [mm] (Seq 2)']);
% xlim(slipRateRange);
% ylim(frictionRange);
% set(gca, 'fontsize', fontsize);
% savename = "../Vitoplots/f_V_x" + num2str(target_x) + "seq2.png";
% print(fig, savename, '-dpng', '-r500');