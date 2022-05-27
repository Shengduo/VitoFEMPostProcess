function res = residualOfStress(x0, Xs, Ys)
% Compute the squared integrated residual stress given X and Y
    res = trapz(Xs, Ys - x0)^2;
end