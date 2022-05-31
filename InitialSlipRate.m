%% Estimate initial slip rate
clc,clear;
close all;

% Parameters
a = 0.011;
b = 0.016;
Drs = 1.5e-6;
Vr = 1e-6;
alpha = 29;
P = 14.3;
si_perturbed = 1.12;
si = P * cosd(alpha) ^ 2 - si_perturbed;
tau = P * cosd(alpha) * sind(alpha);
fr = 0.58;
theta_ini = 0.07;

% Calculate inital slip rate
Vini = exp((tau / si - fr - b * log(Vr * theta_ini / Drs)) / a) * Vr;
disp("Initial slip rate is: " + num2str(Vini) + " [m/s]");