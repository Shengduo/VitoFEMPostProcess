function res = residualOfStress(x0, Xs, Ys, flatflag)
% Compute the squared integrated residual stress given X and Y
    si0 = 14.3 * cosd(29)^2;
    xrange = [-84.9999,  143.6700];
    VS_start = [0.006354, 0.003522, 0.0];
    VS_end = [0.063204, 0.035034, 0];
    Zs = 1e3 * [-0.005, 0.005; -0.0040, 0.0044];
    edge = Zs(1, 2) - Zs(1, 1) - (Zs(2, 2) - Zs(2, 1));
    % parallel_vec = (VS_end - VS_start) ./ norm(VS_end - VS_start, 2);
    faultLength = norm(VS_end - VS_start, 2);
    WirePos1 = [-0.025657, -0.014222, 0];
    interpolate_xs = [xrange(1), norm(VS_start - WirePos1, 2) * 1e3, ...
                  Xs + norm(VS_start - WirePos1, 2) * 1e3, ...
                  (faultLength +  norm(VS_start - WirePos1, 2)) * 1e3, ...
                  xrange(2)];
    interpolate_ys = [x0, x0, Ys, x0, x0];
    if flatflag == 1
        res = (trapz(interpolate_xs, interpolate_ys - x0) * (Zs(2,2) - Zs(2,1))...
               + x0 * (Zs(1, 2) - Zs(1, 1)) * (xrange(2) - xrange(1)) - si0 * (xrange(2) - xrange(1)) * (Zs(1, 2) - Zs(1, 1)))^2;
    else 
        res = (0.5 * trapz(interpolate_xs, interpolate_ys - x0) * (Zs(2,2) - Zs(2,1))...
               + x0 * (Zs(1, 2) - Zs(1, 1)) * (xrange(2) - xrange(1)) - si0 * (xrange(2) - xrange(1)) * (Zs(1, 2) - Zs(1, 1)))^2;
    end
end