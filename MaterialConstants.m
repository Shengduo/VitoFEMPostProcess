% Calculate \lambda, G and density and other elastic constants 
% from cs and cp
clc,clear;
close all;
cs = 1279;
cp = 2662.4;
density = 1200;
la = density * (cp ^ 2 - 2 * cs ^ 2);
G = density * cs ^ 2;
E = G * (3 * la + 2 * G) / (la + G);
nu = la / (2 * (la + G));
s_yy = - 14.3e6;

alpha = 29 / 180 * pi;
tau0 = - s_yy * sin(alpha) * cos(alpha);
sigma0 = - s_yy * cos(alpha) * cos(alpha);

% Apply non-uniform traction
NULoad = 0.75e6;
si_minus = sigma0 - NULoad;
tau_minus = tau0 - NULoad * tan(alpha);

si_plus = sigma0 + NULoad / 3.;
tau_plus = tau0 + NULoad * tan(alpha) / 3.;

ratio = 0.2 * 0.01 / (0.0078 * (0.058832 - 0.006354));

si_plusDiff = sigma0 + NULoad / ratio;
si_minusDiff = sigma0 - NULoad;

tau_plusDiff = tau0 + NULoad * tan(alpha) / ratio;
tau_minusDiff = tau0 - NULoad * tan(alpha);

strain_yy = s_yy / E;
strain_xx = - nu * strain_yy;
strain_zz = - nu * strain_yy;