%% Write the coordinates and corresponding values
clc,clear;
close all;

% Define the ranges for the block, in R^3
ratio = 1.1;
x_range = ratio * [-0.1, 0.1];
y_range = ratio * [-0.1, 0.1];
z_range = ratio * [-0.005, 0.005];

% Other parameters
x_incre = ratio * 0.001;
y_incre = ratio * 0.001;
z_incre = ratio * 0.001;

spaceDim = 3; % Space dimension, default 3
numValues = 3; % [cp, cs, density]

% Normal direction of the plane
normalDir = [-sind(29), cosd(29), 0.];

% Point of the wire position
WirePos1 = [-0.025657; -0.014222; 0]';

% Generate the grid
xs = x_range(1) : x_incre : x_range(2);
ys = y_range(1) : y_incre : y_range(2);
zs = z_range(1) : z_incre : z_range(2);

% Set of points
PointSet = zeros(length(xs) * length(ys) * length(zs), spaceDim);
ValueSet = zeros(length(xs) * length(ys) * length(zs), numValues);
ptId = 1;

% Calculate the point values
for i = 1:1:length(xs)
    for j = 1:1:length(ys)
        for k = 1:1:length(zs)
            PointSet(ptId, :) = [xs(i), ys(j), zs(k)];
            ValueSet(ptId, :) = cal_Values_dist_exp_2dir(PointSet(ptId, :), normalDir, WirePos1);
            ptId = ptId + 1;
        end
    end
end

% Output to a file
dlmwrite('./shitFile.txt', xs, 'delimiter', '\t', 'precision', '%.6f');
dlmwrite('./shitFile.txt', ys, '-append', 'delimiter', '\t', 'precision', '%.6f');
dlmwrite('./shitFile.txt', zs, '-append', 'delimiter', '\t', 'precision', '%.6f');
dlmwrite('./shitFile.txt', [PointSet, ValueSet], '-append', 'Delimiter', '\t', 'precision', '%.6f');
disp(strcat("num-x = ", num2str(length(xs))));
disp(strcat("num-y = ", num2str(length(ys))));
disp(strcat("num-z = ", num2str(length(zs))));
disp(strcat("Number of lines: ", num2str(ptId - 1)));



% Function wrapper cal_Values that gives the values at a distance
function values = cal_Values(pt, normal_dir)
    % Distance is pt dot normal_dir
    dist = abs(dot(pt, normal_dir));
    
    % Compute values at a distance
    values = cal_Values_dist_exp(dist);
end


% Function exponential decay
function values = cal_Values_dist_exp(dist)
    % Dynamic values [cs, cp]
    dyn_values = [1279., 2662.4];
    
    % Quasi-static values [cs, cp]
    static_values = [1080., 1890.];
    
    % Density does not change with space
    rho = 1200.;
    
    % Characteristic length to decay by e
    L = 10.e-3;
    
    % Calculate the values
    values = static_values + (dyn_values - static_values) * exp(- dist / L);
    values = [rho, values];
end

function values = cal_Values_dist_exp_2dir(pt, normal_dir, wire_pos)
    % Normal direction decaying
    L_n = 100.e-3;
    
    % Tangent direction decaying characteristic length
    L_t = 1.e3;

    % Density
    rho = 1200.;

    % Dynamic values [cs, cp]
    dyn_values = [1279., 2662.4];
    % dyn_values = [1279., 2500.];
    
    % Quasi-static values [cs, cp]
    static_values = [1080., 1890.];
    
    % Calculate the values
    dist_n = abs(dot(normal_dir, pt - wire_pos));
    dist_t = sqrt(norm(pt - wire_pos, 2)^2 - dist_n ^ 2);
    values = static_values + (dyn_values - static_values) * exp(- sqrt((dist_n / L_n)^2 + (dist_t / L_t)^2));
    values = [rho, values];
end