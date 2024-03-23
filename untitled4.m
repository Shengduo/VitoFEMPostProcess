% shit
nu = 0.24;
nuu = 0.35;
B = 0.85;
G = 10e9;
c = 4e-8;

alpB = 3*(nuu-nu)/(B*(1+nuu)*(1-2*nu));
kappa = c/(2*G*(1+nu)*B / (3*alpB*(1-alpB*B)*(1-2*nu)));


M = (2 * G * (nuu - nu)) / (alpB ^ 2 * (1 - 2 * nuu) * (1 - 2 * nu));

c = c*((1-nu)*(1-2*nuu) )/( (1-nuu)*(1-2*nu));
M * kappa;

% Have the
V = 1e-7;
theta0 = 3.957060265781263;
DRS = 1.5e-6;

shit = V * theta0 / DRS * log(V * theta0 / DRS);
shit