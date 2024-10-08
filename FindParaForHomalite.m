% Under some constraints, find the parameters for the Homalite
clc,clear;
a1 = 0.011;
b1 = 0.016;
DRS1 = 1.5e-6;
VStar = 1e-6;
fStar = 0.58;
theta1 = 0.006;
fini = tand(29);
V1 = VStar * exp(((fini - fStar) - b1 * log(VStar * theta1 / DRS1)) / a1);

V2 = 1e-7;

% Target vector
targ_vec = [fini, 1e4 * (fStar + a1 * log(1. / VStar) + b1 * log(VStar * theta1 / DRS1)) * DRS1 * log(1. * theta1 / DRS1)];

% Find optimal set of parameters
x0 = [a1, b1, DRS1, theta1];
func = @(x)findPara(x, V2, targ_vec, fStar, VStar);

% Constraints: 0 < a < b, 10e-6 > DRS > 0, theta > 0
A = [-1, 0, 0, 0; ...
     1, -1, 0, 0; ...
     0, 0, -1, 0; ...
     0, 0, 0, -1; ...
     0, 0, 1,  0];

b = [0; 0; -0e-6; 0; 10e-6];


[xmin, fval] = fmincon(func, x0, A, b);
% xmin(1) = xmin(1) - 0.0405;
disp(strcat(["xmin: ", num2str(xmin)]));
disp(strcat(["Minimum value: ", num2str(fval)]));
disp(strcat(["Minimum vector: ", calculate_vec(xmin, V2, fStar ,VStar)]));
disp(strcat(["Minimum vector: ", calculate_vec(xmin, V2, fStar ,VStar)]));

disp(strcat(["Real theta: ", calculate_real_theta(xmin, V2, fStar, VStar, fini)]));

function dist = findPara(x, V, target_vec, fStar, VStar)
    a = x(1);
    b = x(2);
    theta = x(4);
    DRS = x(3);
    vec = [fStar + a * log(V / VStar) + b * log(VStar * theta / DRS), ...
           1e4 * (fStar + a * log(1. / VStar) + b * log(VStar * theta / DRS)) * DRS * log(1. * theta / DRS)];

    dist = norm(vec - target_vec, 2);
end

function vec = calculate_vec(x, V, fStar, VStar)
    a = x(1);
    b = x(2);
    theta = x(4);
    DRS = x(3);
    vec = [fStar + a * log(V / VStar) + b * log(VStar * theta / DRS), ...
           1e4 * (fStar + a * log(1. / VStar) + b * log(VStar * theta / DRS)) * DRS * log(1. * theta / DRS)];
end

function theta_real = calculate_real_theta(x, V, fStar, VStar, f_targ)
    a = x(1);
    b = x(2);
    DRS = x(3);
    theta_real = exp((f_targ - fStar - a * log(V / VStar)) / b) * DRS / VStar;
end