% Calculate the xyz locations for inserting xyz values
setEnvironment;
load('../matFiles/THD58_Friction.mat');
fontsize = 25;

% Find the id for targeted position
target_x = 40;
[val, id] = min(abs(position_THD58_1 - target_x));

% Define ranges for the plots
slipRange = [0, 90];
slipRateRange = [0, 10];
frictionRange = [0, 1];

figNo = 1;
%% f and sliprate vs. slip, seq 1
% Adjust the figure size
fig = figure(figNo);
pos = fig.Position;
pos(4) = pos(4) * 3;
pos(3) = pos(3) * 1.5;
fig.Position = pos;

% f vs. slip
subplot(2, 1, 1);
plot(Slip_avg_THD58_1(id, :), friction_AVG_THD58_1(id, :), 'linewidth', 2.0);
hold on; grid on;
scatter(Slip_avg_THD58_1(id, :), friction_AVG_THD58_1(id, :), 'filled');
xlabel('Slip [$\mu$ m]');
ylabel('Friction coefficient');
title(['Friction vs. slip at x = ', num2str(target_x), ' [mm] (Seq 1)']);
xlim(slipRange);
ylim(frictionRange);
set(gca, 'fontsize', fontsize);

% V vs. slip
subplot(2, 1, 2);
plot(Slip_avg_THD58_1(id, 2:end), Slip_rate_avg_THD58_1(id, :), 'linewidth', 2.0);
hold on; grid on;
scatter(Slip_avg_THD58_1(id, 2:end), Slip_rate_avg_THD58_1(id, :), 'filled');
xlabel('Slip [$\mu$ m]');
ylabel('Slip rate [m/s]');
title(['Slip rate vs. slip at x = ', num2str(target_x), ' [mm] (Seq 1)']);
xlim(slipRange);
ylim(slipRateRange);
set(gca, 'fontsize', fontsize);
savename = "../Vitoplots/fV_slip_x" + num2str(target_x) + "seq1.png";
print(fig, savename, '-dpng', '-r500');


figNo = figNo + 1;
%% f vs slip rate, seq 1
% Adjust the figure size
fig = figure(figNo);
pos = fig.Position;
pos(4) = pos(4) * 1.5;
pos(3) = pos(3) * 1.5;
fig.Position = pos;
plot(Slip_rate_avg_THD58_1(id, :), friction_AVG_THD58_1(id, 2:end), 'linewidth', 2.0);
hold on; grid on;
scatter(Slip_rate_avg_THD58_1(id, :), friction_AVG_THD58_1(id, 2:end), 'filled');
xlabel('Slip rate [m/s]');
ylabel('Friction coefficient');
title(['Friction vs. slip rate at x = ', num2str(target_x), ' [mm] (Seq 1)']);
xlim(slipRateRange);
ylim(frictionRange);
set(gca, 'fontsize', fontsize);
savename = "../Vitoplots/f_V_x" + num2str(target_x) + "seq1.png";
print(fig, savename, '-dpng', '-r500');


figNo = figNo + 1;
%% f and sliprate vs. slip, seq 2
% Adjust the figure size
fig = figure(figNo);
pos = fig.Position;
pos(4) = pos(4) * 3;
pos(3) = pos(3) * 1.5;
fig.Position = pos;

% f vs. slip
subplot(2, 1, 1);
plot(Slip_avg_THD58_2(id, :), friction_AVG_THD58_2(id, :), 'linewidth', 2.0);
hold on; grid on;
scatter(Slip_avg_THD58_2(id, :), friction_AVG_THD58_2(id, :), 'filled');
xlabel('Slip [$\mu$ m]');
ylabel('Friction coefficient');
title(['Friction vs. slip at x = ', num2str(target_x), ' [mm] (Seq 2)']);
xlim(slipRange);
ylim(frictionRange);
set(gca, 'fontsize', fontsize);

% V vs. slip
subplot(2, 1, 2);
plot(Slip_avg_THD58_2(id, 2:end), Slip_rate_avg_THD58_2(id, :), 'linewidth', 2.0);
hold on; grid on;
scatter(Slip_avg_THD58_2(id, 2:end), Slip_rate_avg_THD58_2(id, :), 'filled');
xlabel('Slip [$\mu$ m]');
ylabel('Slip rate [m/s]');
title(['Slip rate vs. slip at x = ', num2str(target_x), ' [mm] (Seq 2)']);
xlim(slipRange);
ylim(slipRateRange);
set(gca, 'fontsize', fontsize);
savename = "../Vitoplots/fV_slip_x" + num2str(target_x) + "seq2.png";
print(fig, savename, '-dpng', '-r500');


figNo = figNo + 1;
%% f vs slip rate, seq 2
% Adjust the figure size
fig = figure(figNo);
pos = fig.Position;
pos(4) = pos(4) * 1.5;
pos(3) = pos(3) * 1.5;
fig.Position = pos;
plot(Slip_rate_avg_THD58_2(id, :), friction_AVG_THD58_2(id, 2:end), 'linewidth', 2.0);
hold on; grid on;
scatter(Slip_rate_avg_THD58_2(id, :), friction_AVG_THD58_2(id, 2:end), 'filled');
xlabel('Slip rate [m/s]');
ylabel('Friction coefficient');
title(['Friction vs. slip rate at x = ', num2str(target_x), ' [mm] (Seq 2)']);
xlim(slipRateRange);
ylim(frictionRange);
set(gca, 'fontsize', fontsize);
savename = "../Vitoplots/f_V_x" + num2str(target_x) + "seq2.png";
print(fig, savename, '-dpng', '-r500');