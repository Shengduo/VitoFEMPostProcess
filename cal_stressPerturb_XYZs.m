% Calculate the xyz locations for inserting xyz values
clc,clear;
close all;
parallel_vec = [0.100000, 0.055430, 0.0] ./ norm([0.100000, 0.055430, 0.0], 2);

% Start of velocity strengthening region
VS_start = [0.006354, 0.003522, 0.0];
WirePos1 = [-0.025657, -0.014222, 0];

% Points to center stress perturbation and length of each perturbed region
% [mm]
points = [8, 18];
interval = 8;
num_points = 11;

% Location of points
XYZs = zeros(size(points, 2) * num_points, 3);
counter = 1;

% Generate XYZs
for pt = 1:1:size(points, 2)
    start = VS_start + (points(pt) - interval / 2) * 1e-3 * parallel_vec;
    for i = 1:1:num_points
        XYZs(counter, :) = start + (i - 1) * (interval / (num_points - 1)) * 1e-3 * parallel_vec;
        counter = counter + 1;
    end
end

%% plot the non-uniform load
figNo = 1;
fig = figure(figNo);
fig.Position(3:4) = 5 * fig.Position(3:4);
NUload = 0;
xrange = [-84.9999,  143.6700];
x_grid = xrange(1) : 1 : xrange(2);
load = 14.3 * cosd(29)^2 * ones(1, size(x_grid, 2));
points = points + norm(VS_start - WirePos1, 2) * 1e3;
XYZloads = zeros(size(XYZs, 1), 3);
for i = 1:1:size(x_grid, 2)
    if (x_grid(i) >= points(1) - interval / 2) && (x_grid(i) <= points(1) + interval / 2)
        load(i) = load(i) - (1 + cos(2 * pi / interval * (x_grid(i) - points(1)))) / 2 * NUload;
    elseif (x_grid(i) >= points(2) - interval / 2) && (x_grid(i) <= points(2) + interval / 2)
        load(i) = load(i) + (1 + cos(2 * pi / interval * (x_grid(i) - points(2)))) / 2 * NUload;
    end
end

XYZnorms = zeros(1, size(XYZs, 1));
points = [8, 18];
points = points + norm(VS_start, 2) * 1e3;
for i = 1:1:size(XYZloads, 1)
    XYZnorms(i) = 1e3 * norm(XYZs(i, :), 2);
    if (XYZnorms(i) >= points(1) - interval / 2) && (XYZnorms(i) <= points(1) + interval / 2)
        XYZloads(i, 3) = (1 + cos(2 * pi / interval * (XYZnorms(i) - points(1)))) / 2 * NUload;
    elseif (XYZnorms(i) >= points(2) - interval / 2) && (XYZnorms(i) <= points(2) + interval / 2)
        XYZloads(i, 3) = -(1 + cos(2 * pi / interval * (XYZnorms(i) - points(2)))) / 2 * NUload;
    end
end

plot(x_grid, load, 'linewidth', 2.0);
hold on; grid on;
xlabel('Distance along the fault [mm]', 'interpreter', 'latex');
ylabel({'Initial normal', 'stress [MPa]'}, 'Interpreter', 'latex');
title('Distribution of initial normal stress along the fault');
axis equal;
xlim(xrange);
ylim([5, 17]);
set(gca, 'fontsize', 25);

% Write into files
% Write changable parameters into a '.txt' file
txtname = "XYZs.txt";
fileID = fopen(txtname, 'w');

Zs = [-0.0051, 0.0051; -0.0041, 0.0045; -0.0039, 0.0043];


% Write into the files
for i = 1:1:size(XYZs, 1)
    % Write 9 lines
    for shit = 1:1:3
        fprintf(fileID, '%9s', num2str(XYZs(i, 1), '%6f')); 
        fprintf(fileID, '%10s', num2str(0.1, '%6f')); 
        fprintf(fileID, '%10s', num2str(Zs(shit, 1), '%6f'));
        if shit == 3
            fprintf(fileID, '%10s', num2str(XYZloads(i, 1), '%6f')); 
            fprintf(fileID, '%10s', num2str(XYZloads(i, 2), '%6f')); 
            fprintf(fileID, '%10s', num2str(XYZloads(i, 3), '%6f'));
            fprintf(fileID, '\n');
        else
            fprintf(fileID, '%10s', num2str(0.0, '%6f')); 
            fprintf(fileID, '%10s', num2str(0.0, '%6f')); 
            fprintf(fileID, '%10s', num2str(0.0, '%6f'));
            fprintf(fileID, '\n');
        end
        
        % ================================================
        fprintf(fileID, '%9s', num2str(XYZs(i, 1), '%6f')); 
        fprintf(fileID, '%10s', num2str(0.1, '%6f')); 
        fprintf(fileID, '%10s', num2str(Zs(shit, 2), '%6f'));
        if shit == 3
            fprintf(fileID, '%10s', num2str(XYZloads(i, 1), '%6f')); 
            fprintf(fileID, '%10s', num2str(XYZloads(i, 2), '%6f')); 
            fprintf(fileID, '%10s', num2str(XYZloads(i, 3), '%6f'));
            fprintf(fileID, '\n');
        else
            fprintf(fileID, '%10s', num2str(0.0, '%6f')); 
            fprintf(fileID, '%10s', num2str(0.0, '%6f')); 
            fprintf(fileID, '%10s', num2str(0.0, '%6f'));
            fprintf(fileID, '\n');
        end
        
        % ================================================
        fprintf(fileID, '%9s', num2str(XYZs(i, 1), '%6f')); 
        fprintf(fileID, '%10s', num2str(-0.1, '%6f')); 
        fprintf(fileID, '%10s', num2str(Zs(shit, 1), '%6f'));
        if shit == 3
            fprintf(fileID, '%10s', num2str(XYZloads(i, 1), '%6f')); 
            fprintf(fileID, '%10s', num2str(XYZloads(i, 2), '%6f')); 
            fprintf(fileID, '%10s', num2str(XYZloads(i, 3), '%6f'));
            fprintf(fileID, '\n');
        else
            fprintf(fileID, '%10s', num2str(0.0, '%6f')); 
            fprintf(fileID, '%10s', num2str(0.0, '%6f')); 
            fprintf(fileID, '%10s', num2str(0.0, '%6f'));
            fprintf(fileID, '\n');
        end
        
        % ================================================
        fprintf(fileID, '%9s', num2str(XYZs(i, 1), '%6f')); 
        fprintf(fileID, '%10s', num2str(-0.1, '%6f')); 
        fprintf(fileID, '%10s', num2str(Zs(shit, 2), '%6f'));
        if shit == 3
            fprintf(fileID, '%10s', num2str(XYZloads(i, 1), '%6f')); 
            fprintf(fileID, '%10s', num2str(XYZloads(i, 2), '%6f')); 
            fprintf(fileID, '%10s', num2str(XYZloads(i, 3), '%6f'));
            fprintf(fileID, '\n');
        else
            fprintf(fileID, '%10s', num2str(0.0, '%6f')); 
            fprintf(fileID, '%10s', num2str(0.0, '%6f')); 
            fprintf(fileID, '%10s', num2str(0.0, '%6f'));
            fprintf(fileID, '\n');
        end
        fprintf(fileID, '\n');
    end
end

% Output XYZ
for i = 1:1:size(XYZs, 1)
    fprintf(fileID, '%10s', num2str(XYZs(i, 1), '%6f')); 
end
fclose(fileID);
