clc,clear; 
close all;
Vflag = 1;

% Get peak friction coefficient from files
load('../matFiles/THD58_Friction.mat');
peakFs = zeros(2, size(friction_AVG_THD58_1, 1));
peakVs = zeros(2, size(friction_AVG_THD58_1, 1));

% File 1
[F, ind] = max(friction_AVG_THD58_1, [], 2);
peakFs(1, :) = F';
for i = 1:1:size(Slip_rate_avg_THD58_1, 1)
    peakVs(1, i) = Slip_rate_avg_THD58_1(i, min(ind(i),size(Slip_rate_avg_THD58_1, 2)));
end

% File 2
[F, ind] = max(friction_AVG_THD58_2, [], 2);
peakFs(2, :) = F';
for i = 1:1:size(Slip_rate_avg_THD58_2, 1)
    peakVs(2, i) = Slip_rate_avg_THD58_2(i, max(ind(i),size(Slip_rate_avg_THD58_2, 2)));
end
Xs = position_THD58_2;

% Calculate change of delta_f
P = 14.3;
alpha = 29;
f_0 = tand(alpha);
si0 = P * cosd(alpha)^2;

deltaFs = peakFs - f_0;
Vini = 1e-9;
Vdyn = 1;

a = 0.011;

% Compute real normal stress around
si = deltaFs * si0 / (a * log(Vdyn / Vini));
if Vflag == 1
    si = deltaFs * si0 ./ (a .* log(abs(peakVs) ./ Vini));
end


% Plot distribution of initial normal stress
fig = figure(1);
plot(Xs, si(1, :), '--', 'color', '#0072BD', 'linewidth', 3.0);
hold on; grid on;
plot(Xs, si(2, :), '--', 'color', '#D95319', 'linewidth', 3.0);

si_smooth = si;
span = 0.3;
si_smooth(1, :) = smooth(si(1, :), span, 'sgolay');
si_smooth(2, :) = smooth(si(2, :), span, 'sgolay');

plot(Xs, si_smooth(1, :), 'color', '#0072BD', 'linewidth', 2.0);
hold on; grid on;
plot(Xs, si_smooth(2, :), 'color', '#D95319', 'linewidth', 2.0);

xlim([0, 50]); ylim([0, 25]);
xlabel('Distance Along The Fault [mm]', 'interpreter', 'latex');
ylabel('Actual Normal Stress [MPa]', 'interpreter', 'latex');
legend('1st explosion', '2nd explosion', '1st fit', '2nd fit', 'location', 'best', 'interpreter', 'latex');
set(gca, 'fontsize', 15);

%Save the figure
print(fig ,'../Vitoplots/sigmaDistri.png', '-dpng', '-r500');

% Save the data
save('../matFiles/realNormalStress.mat', 'Xs', 'si', 'si_smooth');

%% Plot seperately for sequence 1 and 2
fig = figure(2);
plot(Xs, si(1, :), '--', 'color', '#0072BD', 'linewidth', 3.0);
hold on; grid on;
plot(Xs, si_smooth(1, :), 'color', '#0072BD', 'linewidth', 2.0);
xlim([0, 50]); ylim([0, 25]);
yline(si0, '-.k', 'linewidth', 1.5);
xlabel('Distance Along The Fault [mm]', 'interpreter', 'latex');
ylabel('Actual Normal Stress [MPa]', 'interpreter', 'latex');
legend('1st explosion', '1st fit', 'Uniform $\sigma_0$', 'location', 'best', 'interpreter', 'latex');
set(gca, 'fontsize', 15);
print(fig ,'../Vitoplots/sigmaDistri-1.png', '-dpng', '-r500');

fig = figure(3);
plot(Xs, si(2, :), '--', 'color', '#D95319', 'linewidth', 3.0);
hold on; grid on;
plot(Xs, si_smooth(2, :), 'color', '#D95319', 'linewidth', 2.0);
xlim([0, 50]); ylim([0, 25]);
yline(si0, '-.k', 'linewidth', 1.5);
xlabel('Distance Along The Fault [mm]', 'interpreter', 'latex');
ylabel('Actual Normal Stress [MPa]', 'interpreter', 'latex');
legend('2nd explosion', '2nd fit', 'Uniform $\sigma_0$', 'location', 'best', 'interpreter', 'latex');
set(gca, 'fontsize', 15);
print(fig ,'../Vitoplots/sigmaDistri-2.png', '-dpng', '-r500');