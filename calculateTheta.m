clc,clear; 
close all;

Vini = 1e-6;
f = tand(29);
A = 0.011;
B = 0.006;
DRS = 1.5e-6;
Vprev9 = Vini * exp(- 0.006 / 0.011 * 6 * log(10.));

temp9 = f - 0.58 - A * log(Vprev9 / 1e-6);



theta2 = exp(temp9 / B) * DRS * 1e6;

Vprev6 = Vprev9 * 1e3;
tempPrev = f - 0.58 - A * log(Vprev6 / 1e-6);
theta2_prev = exp(tempPrev / B) * DRS * 1e6;

bnew = log(1e-6 * theta2_prev / 1.5e-6) / log(1e-6 * theta2 / 1.5e-6) * 0.006


% ---------------------------------------------------------------------------
f = tand(29);
A = 0.003;
B = 0.008;
DRS = 1.5e-6;
fStar = 0.58;
VStar = 1e-6;
V_targ = 1e-7;

theta_ini = exp((f - fStar - A * log(V_targ / VStar)) / B) * DRS / VStar
