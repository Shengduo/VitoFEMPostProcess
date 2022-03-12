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

strain_yy = s_yy / E;
strain_xx = - nu * strain_yy;
strain_zz = - nu * strain_yy;