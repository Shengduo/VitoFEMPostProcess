function [Gamma, si_xy] = FracSigmaGivenX(x, y, cf, cs, cd, rho, KII)
    %% Following Appendix A of Kammer and McLaskey (2019)
    % Fracture estimates from large-scale laboratory earthquakes

    % x - horizontal distance from the crack tip
    % y - vertical distance above the crack path
    % cf - crack propagating speed
    % KII - Stress intensity factor
    % Return si_xy - Computed [si_xy]
    %        Gamma - Energy release rate
    
    % Geometry
    theta = atan(y ./ x);
    theta(theta < 0) = theta(theta < 0) + pi;

    r = sqrt(x .^ 2 + y .^ 2);

    alp_s = sqrt(1 - (cf ./ cs) .^ 2);
    alp_d = sqrt(1 - (cf ./ cd) .^ 2);
    
    % Rayleigh function
    D_cf = 4 .* alp_s .* alp_d - (1 + alp_s .^ 2) .^ 2;

    
    theta_s = atan(alp_s .* tan(theta));
    theta_d = atan(alp_s .* tan(theta));
    
    theta_s(theta_s < 0) = theta_s(theta_s < 0) + pi;
    theta_d(theta_d < 0) = theta_d(theta_d < 0) + pi;
    

    gamma_s = sqrt(1 - cf / cs .* sin(theta));
    gamma_d = sqrt(1 - cf / cd .* sin(theta));
    
    % Compute si_xy
    Sigma_xy = (4 * alp_s * alp_d .* cos(theta_d / 2) ./ sqrt(gamma_d) - ...
                (1 + alp_s ^ 2) ^ 2 .* cos(theta_s / 2) ./ sqrt(gamma_s)) ./ D_cf;
    si_xy = KII ./ sqrt(2 * pi .* r) .* Sigma_xy;
    
    % Compute EnergyReleaseRate
    Gamma = computeEnergyReleaseRate(cf, cs, cd, rho, KII);
end


%% Function that computes the energy relaase rate
function Gamma = computeEnergyReleaseRate(cf, cs, cd, rho, KII)
    % Elastic constants
    la = rho * (cd ^ 2 - 2 * cs ^ 2);
    G = rho * cs ^ 2;
    E = G * (3 * la + 2 * G) / (la + G);
    nu = la / (2 * (la + G));

    % Geometry

    alp_s = sqrt(1 - (cf ./ cs) .^ 2);
    alp_d = sqrt(1 - (cf ./ cd) .^ 2);
    
    % Rayleigh function
    D_cf = 4 .* alp_s .* alp_d - (1 + alp_s .^ 2) .^ 2;
    
    % Compute energy release rate
    fIICf = alp_s / (1 - nu) / D_cf * (cf / cs) ^ 2;
    Gamma = (1 - nu ^ 2) / E * KII ^ 2 * fIICf;
end

