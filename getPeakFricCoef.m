clc,clear; close all;
% Get peak friction coefficient from files
files = ["../matFiles/THD58_1_Friction.mat", "../matFiles/THD58_2_Friction.mat"];

for i = 1:1:size(files, 2)
    load(files(i));
    peakFs(i, :) = peakFricCoef';
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

% Plot distribution of initial normal stress
fig = figure(1);
plot(Xs, si(1, :), 'linewidth', 2.0);
hold on; grid on;
plot(Xs, si(2, :), 'linewidth', 2.0);
legend('First explosion', 'Second explosion', 'location', 'best');
set(gca, 'fontsize', 20);

%Save the figure
print(fig ,'../matFiles/sigmaDistri.png', '-dpng', '-r500');

% Save the data
save('../matFiles/realNormalStress.mat', 'Xs', 'si');