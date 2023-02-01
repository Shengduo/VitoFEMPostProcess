% Calculate the xyz locations for inserting xyz values
setEnvironment;
load('../matFiles/THD58_Friction.mat');
load('../matFiles/realNormalStress.mat');
fontsize = 25;

% Define ranges for the plots
slipRange = [0, 40];
slipRateRange = [0, 10];
frictionRange = [0, 1];

% Compute si0, tau0 [MPa]
P = 14.3;
si0 = P * cosd(29)^2;
tau0 = P * cosd(29) * sind(29);

% Compute change of shear stress
shearChange_seq1 = friction_AVG_THD58_1 * si0 - tau0;
% shearChangeRange = 1.2 * [min(min(shearChange_seq1)), max(max(shearChange_seq1))];
shearChangeRange = [-4, 4];

% Set the plotting positions
locs = 4:3:30;

% Plot and save tauChange vs. sliprate and slip
for i = 1:1:length(locs)
    % Find the index of location
    [~, idx] = min(abs(position_THD58_1 - locs(i))); 
    disp(idx);
%     figNo = 1;
%     %% f and sliprate vs. slip, seq 1
%     % Adjust the figure size
%     fig = figure(figNo);
%     pos = fig.Position;
%     pos(4) = pos(4) * 3.0;
%     pos(3) = pos(3) * 1.5;
%     fig.Position = pos;

    % delta tau vs. slip
    subplot(1, length(locs), i);
    plot(Slip_avg_THD58_1(idx, :), shearChange_seq1(idx, :), 'linewidth', 2.0);
    hold on; grid on;
    scatter(Slip_avg_THD58_1(idx, :), shearChange_seq1(idx, :), 'filled');
    xlabel('Slip $[\mathrm{\mu m}]$');
    if i == 1
        ylabel('$\Delta \tau$ [MPa]');
    end
    title(strcat('$x_1 = $', num2str(position_THD58_1(idx))));

    xlim(slipRange);
    ylim(shearChangeRange);
    set(gca, 'fontsize', fontsize);

%     % V vs. slip
%     subplot(2, 1, 2);
%     plot(Slip_avg_THD58_1(id, 2:end), Slip_rate_avg_THD58_1(id, :), 'linewidth', 2.0);
%     hold on; grid on;
%     scatter(Slip_avg_THD58_1(id, 2:end), Slip_rate_avg_THD58_1(id, :), 'filled');
%     xlabel('Slip [$\mu$ m]');
%     ylabel('Slip rate [m/s]');
%     title(['Slip rate vs. slip at x = ', num2str(target_x), ' [mm] (Seq 1)']);
%     xlim(slipRange);
%     ylim(slipRateRange);
%     set(gca, 'fontsize', fontsize);
%     savename = "../plots20221117/experiment" + num2str() + "seq1.png";
%     print(fig, savename, '-dpng', '-r500');


%     figNo = figNo + 1;
%     %% f vs slip rate, seq 1
%     % Adjust the figure size
%     fig = figure(figNo);
%     pos = fig.Position;
%     pos(4) = pos(4) * 1.5;
%     pos(3) = pos(3) * 1.5;
%     fig.Position = pos;
%     plot(Slip_rate_avg_THD58_1(id, :), friction_AVG_THD58_1(id, 2:end), 'linewidth', 2.0);
%     hold on; grid on;
%     scatter(Slip_rate_avg_THD58_1(id, :), friction_AVG_THD58_1(id, 2:end), 'filled');
%     xlabel('Slip rate [m/s]');
%     ylabel('Friction coefficient');
%     title(['Friction vs. slip rate at x = ', num2str(target_x), ' [mm] (Seq 1)']);
%     xlim(slipRateRange);
%     ylim(frictionRange);
%     set(gca, 'fontsize', fontsize);
%     savename = "../Vitoplots/f_V_x" + num2str(target_x) + "seq1.png";
%     print(fig, savename, '-dpng', '-r500');
% 
% 
%     figNo = figNo + 1;
% 
%     close all;
end
sgtitle('Experiment $\Delta \tau$ vs. Slip', 'fontsize', 25);

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