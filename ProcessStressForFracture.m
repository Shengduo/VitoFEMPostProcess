clc,clear;
close all;
setEnvironment;

% Load the file
mat_name = "../Stress-X-T_mats/W8_4ExcCorr2NPDirWithWallDRS1.5_1.5ModA0.016AmB0.014Load5_Vw2_fw0.1_theta0.1434_-7_NULoad2dir0_duration200_8.mat";
load(mat_name);

%% Calculate the left branch cf speed
WireX = -norm([0.006354, 0.003522, 0]' - [-0.025657; -0.014222; 0], 2);
leftBr_idx = (X_save <= WireX);
rightBr_idx = (X_save > WireX);

% Turn on for left branch only
V = V_save(leftBr_idx, :);
X = X_save(leftBr_idx);
si_xy = si_xy_save(leftBr_idx, :);

%% Plot sliprate-x profile
plotFlag = true;
if plotFlag == true
    figure(1);
    for i = 20:2:40
        plot(X, V(:, i) + 1e6 * T_save(i));
        hold on;
    end
    xlabel("$x_1$", "FontSize", 20);
    ylabel("Time [$\mathrm{\mu s}$]", "FontSize", 20);
    title("Evolution of $V$", "FontSize", 20)
    set(gca, "FontSize", 20);
end

%% Plot delta_si_xy-x profile
plotFlag = true;
if plotFlag == true
    figure(2);
    for i = 20:2:40
        plot(X, si_xy(:, i) / 1e6 + 1e6 * T_save(i));
        hold on;
    end
    xlabel("$x_1$", "FontSize", 20);
    ylabel("Time [$\mathrm{\mu s}$]", "FontSize", 20);
    title("Evolution of $\Delta \sigma_{xy}$", "FontSize", 20)
    set(gca, "FontSize", 20);
end


%% By maximum velocity criteria
xLength = 10e-3;
[Vmax, IVmax] = max(V);
xTips = X(IVmax);
for i = 20:1:40
    xMinusXTip = X - xTips(i);
    % idx = ((xMinusXTip > -xLength / 2) & (xMinusXTip < xLength / 2));
    idx = IVmax(i) - 150 : 1 : IVmax(i) + 50;
    temp = xMinusXTip(idx);
    tempV = V(idx, i);
    tempSi = si_xy(idx, i);
    if i >= 41
        temp = temp(1:size(xMinusXTips, 2));
        tempV = tempV(1:size(xMinusXTips, 2));
        tempSi = tempSi(1:size(xMinusXTips, 2));
    end
    xMinusXTips(i - 19, :) = temp;
    V_aligned(i - 19, :) = tempV;
    si_xy_aligned(i - 19, :) = tempSi;
end

%% Plot aligned sliprate-x-x_tip profile
plotFlag = true;
if plotFlag == true
    figure(3);
    for i = 1:1:21
        plot(xMinusXTips(i, :), V_aligned(i, :) + 1e6 * T_save(i + 19));
        hold on;
    end
    xlabel("$x_1 - x_{tip}$", "FontSize", 20);
    ylabel("Time [$\mathrm{\mu s}$]", "FontSize", 20);
    title("Evolution of $V$", "FontSize", 20)
    set(gca, "FontSize", 20);
end

%% Plot aligned delta_si_xy-x_tip profile
plotFlag = true;
if plotFlag == true
    figure(4);
    for i = 1:1:21
        % plot(xMinusXTips(i, :), si_xy_aligned(i, :) / 1e6 + 1e6 * T_save(i + 19));
        plot(xMinusXTips(i, :), si_xy_aligned(i, :) / 1e6);
        hold on;
    end
    xlabel("$x_1 - x_{tip}$", "FontSize", 20);
    % ylabel("Time [$\mathrm{\mu s}$]", "FontSize", 20);
    ylabel("$\Delta \sigma_{xy}$ [$\mathrm{\mu s}$]", "FontSize", 20);
    ylim([0., 2.]);
    title("Evolution of $\Delta \sigma_{xy}$", "FontSize", 20)
    set(gca, "FontSize", 20);
end


%% Calculate wave speed
cf = zeros(1, size(T_save, 2));
for i = 11:1:size(T_save, 2) - 10
    cf(i) = abs((xTips(i + 10) - xTips(i - 10)) / (T_save(i + 10) - T_save(i - 10)));
end

cf_Interested = cf(20:1:40);

%% Fit for Gamma
cf_Interested_mean = mean(cf_Interested);
xMinusXTip_mean = -mean(xMinusXTips, 1);
y = 3.5e-3;
si_xy_aligned_mean = mean(si_xy_aligned, 1);
V_aligned_mean = mean(V_aligned, 1);

cd = 2662.4;
cs = 1279;
rho = 1200;

% Find KII
KII0 = 1.;
fun = @(KII) errAtKII(xMinusXTip_mean, y, cf_Interested_mean, cs, cd, rho, KII, si_xy_aligned_mean);
options = optimset('TolFun', 1e-12, 'TolX', 1e-2);
[KII_fitted, fval] = fminsearch(fun, KII0, options);
[Gamma_fitted, si_xy_fitted] = FracSigmaGivenX(xMinusXTip_mean, y, cf_Interested_mean, cs, cd, rho, KII_fitted);
shit = fun(1.1 * KII_fitted);
%% Plot fitted curve of si
plotFlag = true;
if plotFlag == true
    figure(5);
    plot(xMinusXTip_mean, si_xy_aligned_mean / 1e6, 'LineWidth', 3.0);
    hold on;
    plot(xMinusXTip_mean, si_xy_fitted / 1e6, 'LineWidth', 1.5);
    % plot(xMinusXTip_mean, 1.5 * si_xy_fitted / 1e6, 'LineWidth', 1.5);
    legend("Simulation", "LEFM-fitted", 'Location', 'best');
    title(strcat("$\Gamma = $", num2str(Gamma_fitted), "$\ \mathrm{[J/m^2]}$"));
    xlabel("$x-x_{tip}$ [m]");
    ylabel("$\Delta \sigma_{xy}$ [MPa]");
    ylim([-1.5, 1.5]);
    set(gca, 'Fontsize', 20);
    set(gcf,'color','w');
end

% [Gamma, si_xy] = FracSigmaGivenX(x, y, cf, cs, cd, rho, KII);

function err = errAtKII(x, y, cf, cs, cd, rho, KII, si_xy_targ)
    [~, si_xy] = FracSigmaGivenX(x, y, cf, cs, cd, rho, KII);
    err = dot(si_xy_targ - si_xy, si_xy_targ - si_xy);
end
