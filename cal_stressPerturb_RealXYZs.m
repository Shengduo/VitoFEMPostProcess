% Calculate the xyz locations for inserting xyz values
clc,clear;
close all;
load('../matFiles/realNormalStress.mat');

% Which sequence to use
seq_ID = 1;
yToxRatio = 2;
flatflag = 1;

% Start of velocity strengthening region
VS_start = [0.006354, 0.003522, 0.0];
VS_end = [0.063204, 0.035034, 0];
parallel_vec = (VS_end - VS_start) ./ norm(VS_end - VS_start, 2);
faultLength = norm(VS_end - VS_start, 2);
WirePos1 = [-0.025657, -0.014222, 0];

% Points to center stress perturbation and length of each perturbed region
% [mm]
num_points = 50;

% Location of points
XYZs = zeros(num_points, 3);
counter = 1;

% Generate XYZs
usableFault = [min(Xs), max(Xs)];
usableFaultLength = usableFault(2) - usableFault(1);
XYZs1D = zeros(1, num_points);
for pt = 1:1:num_points
    XYZs(pt, :) = VS_start + ((pt - 1) * usableFaultLength / (num_points - 1) + usableFault(1)) * parallel_vec* 1e-3;
    XYZs1D(1, pt) = ((pt - 1) * usableFaultLength / (num_points - 1) + usableFault(1)) + norm(VS_start - WirePos1, 2) * 1e3;
end
XYZs = [VS_start; XYZs; VS_end];
XYZs1D = [norm(VS_start - WirePos1, 2) * 1e3, XYZs1D, norm(VS_end - WirePos1, 2) * 1e3];

%% plot the non-uniform load
xrange = [-84.9999,  143.6700];
si0 = 14.3 * cosd(29)^2;
interpolate_xs = [xrange(1), norm(VS_start - WirePos1, 2) * 1e3, ...
                  Xs + norm(VS_start - WirePos1, 2) * 1e3, ...
                  (faultLength +  norm(VS_start - WirePos1, 2)) * 1e3, ...
                  xrange(2)];
Xs_wire = Xs + norm(VS_start - WirePos1, 2) * 1e3;

% Find the background stress
[si0_1, Fval1] = fminsearch(@(x0) residualOfStress(x0, Xs, si_smooth(1, :), flatflag), si0);
[si0_2, Fval2] = fminsearch(@(x0) residualOfStress(x0, Xs, si_smooth(2, :), flatflag), si0);
disp(["Background normal stress of sequence 1 is: ", num2str(si0_1), " MPa"])
disp(["Residual stress of sequence 1 is: ", num2str(Fval1), " MPa"]);
disp(["Background normal stress of sequence 2 is: ", num2str(si0_2), " MPa"])
disp(["Residual stress of sequence 2 is: ", num2str(Fval2), " MPa"]);

interpolate_ys(1, :) = [si0_1, si0_1, si_smooth(1, :), si0_1, si0_1];
interpolate_ys(2, :) = [si0_2, si0_2, si_smooth(2, :), si0_2, si0_2];

x_grid = xrange(1) : 1 : xrange(2);
Load(1, :) = interp1(interpolate_xs, interpolate_ys(1, :), x_grid);
Load(2, :) = interp1(interpolate_xs, interpolate_ys(2, :), x_grid);

XYZloads = zeros(size(XYZs, 1), 3);
xrange = xrange - 1e3 * norm(VS_start - WirePos1, 2);

figNo = 1;
%% Plot for sequence 1
fig = figure(figNo);
fig.Position(4:4) = 2 * fig.Position(4:4);
fig.Position(3:3) = 4 * fig.Position(3:3);
plot(x_grid - 1e3 * norm(VS_start - WirePos1, 2), Load(1, :), 'linewidth', 2.0);
hold on; grid on;
yline(si0, '-.k', 'linewidth', 1.5);
xline(norm(VS_start - VS_start, 2) * 1e3, 'r', 'linewidth', 1.5);
xline(norm(VS_end - VS_start, 2) * 1e3, 'r', 'linewidth', 1.5);
text(norm(VS_start - VS_start, 2) * 1e3 + 20, 3, 'Region 2', 'color', 'r', 'Fontsize', 20);

xlabel('Distance along the fault [mm]', 'interpreter', 'latex');
ylabel({'Initial normal', 'stress [MPa]'}, 'Interpreter', 'latex');
title('Distribution of initial normal stress along the fault infered by sequence 1');
daspect([yToxRatio, 1, 1]);
xlim(xrange);
ylim([0, 20]);
set(gca, 'fontsize', 25);
if flatflag == 1
    print(fig ,'../PaperPlots/sigmaDistriWholeFault-1.png', '-dpng', '-r500');
else
    print(fig ,'../PaperPlots/sigmaDistriWholeFault-1_nonflat.png', '-dpng', '-r500');
end
figNo = figNo + 1;

%% Plot for sequence 2
fig = figure(figNo);
fig.Position(4:4) = 2 * fig.Position(4:4);
fig.Position(3:3) = 4 * fig.Position(3:3);
plot(x_grid - 1e3 * norm(VS_start - WirePos1, 2), Load(2, :), 'color', '#D95319', 'linewidth', 2.0);
hold on; grid on;
yline(si0, '-.k', 'linewidth', 1.5);
xline(norm(VS_start - VS_start, 2) * 1e3, 'r', 'linewidth', 1.5);
xline(norm(VS_end - VS_start, 2) * 1e3, 'r', 'linewidth', 1.5);
text(norm(VS_start - VS_start, 2) * 1e3 + 20, 3, 'Region 2', 'color', 'r', 'Fontsize', 20);

xlabel('Distance along the fault [mm]', 'interpreter', 'latex');
ylabel({'Initial normal', 'stress [MPa]'}, 'Interpreter', 'latex');
title('Distribution of initial normal stress along the fault infered by sequence 2');
daspect([yToxRatio, 1, 1]);
xlim(xrange);
ylim([0, 20]);
set(gca, 'fontsize', 25);
if flatflag == 1
    print(fig ,'../PaperPlots/sigmaDistriWholeFault-2.png', '-dpng', '-r500');
else
    print(fig ,'../PaperPlots/sigmaDistriWholeFault-2_nonflat.png', '-dpng', '-r500');
end
figNo = figNo + 1;


%% Write into files
% Compute the stress at certain points
XYZloads = zeros(size(XYZs, 1), 3);
XYZloads(:, 3) = si0 - interp1(interpolate_xs, interpolate_ys(seq_ID, :), XYZs1D)';

% Add first and last XYZ
XYZs = [-0.110000, 0.0, 0.0; XYZs; 0.110000, 0.0, 0.0];
XYZloads = [XYZloads(1, :); XYZloads; XYZloads(end, :)];

% Write changable parameters into a '.txt' file
if flatflag == 1
    %% Flat file
    txtname = "../matFiles/XYZs.txt";
    fileID = fopen(txtname, 'w');
    Zs = [-0.0051, 0.0051; -0.0041, 0.0045; -0.0039, 0.0043];
    si0_reg = si0 - [si0_1, si0_2];

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
                fprintf(fileID, '%10s', num2str(si0_reg(seq_ID), '%6f'));
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
                fprintf(fileID, '%10s', num2str(si0_reg(seq_ID), '%6f'));
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
                fprintf(fileID, '%10s', num2str(si0_reg(seq_ID), '%6f'));
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
                fprintf(fileID, '%10s', num2str(si0_reg(seq_ID), '%6f'));
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
else
    %% Non-flat
    txtname = "../matFiles/XYZs_nonflat.txt";
    fileID = fopen(txtname, 'w');
    Zs = [-0.0051, 0.0051; -0.0041, 0.0045; -0.0039, 0.0043];
    Zmid = 0.0002;
    si0_reg = si0 - [si0_1, si0_2];

    % Write into the files
    for i = 1:1:size(XYZs, 1)
        % Write 9 lines
        for shit = 1:1:3
            fprintf(fileID, '%9s', num2str(XYZs(i, 1), '%6f')); 
            fprintf(fileID, '%10s', num2str(0.1, '%6f')); 
            fprintf(fileID, '%10s', num2str(Zs(shit, 1), '%6f'));
            fprintf(fileID, '%10s', num2str(0.0, '%6f')); 
            fprintf(fileID, '%10s', num2str(0.0, '%6f')); 
            fprintf(fileID, '%10s', num2str(si0_reg(seq_ID), '%6f'));
            fprintf(fileID, '\n');

            % ================================================
            fprintf(fileID, '%9s', num2str(XYZs(i, 1), '%6f')); 
            fprintf(fileID, '%10s', num2str(0.1, '%6f')); 
            fprintf(fileID, '%10s', num2str(Zs(shit, 2), '%6f'));
            fprintf(fileID, '%10s', num2str(0.0, '%6f')); 
            fprintf(fileID, '%10s', num2str(0.0, '%6f')); 
            fprintf(fileID, '%10s', num2str(si0_reg(seq_ID), '%6f'));
            fprintf(fileID, '\n');

            % ================================================
            fprintf(fileID, '%9s', num2str(XYZs(i, 1), '%6f')); 
            fprintf(fileID, '%10s', num2str(-0.1, '%6f')); 
            fprintf(fileID, '%10s', num2str(Zs(shit, 1), '%6f'));
            fprintf(fileID, '%10s', num2str(0.0, '%6f')); 
            fprintf(fileID, '%10s', num2str(0.0, '%6f')); 
            fprintf(fileID, '%10s', num2str(si0_reg(seq_ID), '%6f'));
            fprintf(fileID, '\n');

            % ================================================
            fprintf(fileID, '%9s', num2str(XYZs(i, 1), '%6f')); 
            fprintf(fileID, '%10s', num2str(-0.1, '%6f')); 
            fprintf(fileID, '%10s', num2str(Zs(shit, 2), '%6f'));
            fprintf(fileID, '%10s', num2str(0.0, '%6f')); 
            fprintf(fileID, '%10s', num2str(0.0, '%6f')); 
            fprintf(fileID, '%10s', num2str(si0_reg(seq_ID), '%6f'));
            fprintf(fileID, '\n');
            fprintf(fileID, '\n');
        end
        % ====================================================
        fprintf(fileID, '%9s', num2str(XYZs(i, 1), '%6f')); 
        fprintf(fileID, '%10s', num2str(0.1, '%6f')); 
        fprintf(fileID, '%10s', num2str(Zmid, '%6f'));
        fprintf(fileID, '%10s', num2str(0.0, '%6f')); 
        fprintf(fileID, '%10s', num2str(0.0, '%6f')); 
        fprintf(fileID, '%10s', num2str(XYZloads(i, 3), '%6f'));
        fprintf(fileID, '\n');

        % ================================================
        fprintf(fileID, '%9s', num2str(XYZs(i, 1), '%6f')); 
        fprintf(fileID, '%10s', num2str(-0.1, '%6f')); 
        fprintf(fileID, '%10s', num2str(Zmid, '%6f'));
        fprintf(fileID, '%10s', num2str(0.0, '%6f')); 
        fprintf(fileID, '%10s', num2str(0.0, '%6f')); 
        fprintf(fileID, '%10s', num2str(XYZloads(i, 3), '%6f'));
        fprintf(fileID, '\n');
        fprintf(fileID, '\n');
    end

    % Output XYZ
    for i = 1:1:size(XYZs, 1)
        fprintf(fileID, '%10s', num2str(XYZs(i, 1), '%6f')); 
    end
    fclose(fileID);
end


