function [levels, t, scaledErr] = sch_1d_exact_err(tmax, lmin, lmax, lambda, idtype, idpar, vtype, vpar)
% Calculates the error between the 1D Schrodinger Eq. and the exact
% solution
% 
% Inputs
% tmax:      Maximum integration time
% level:     Discretization level
% lambda:    dt/dx
% idtype:    Selects initial condition type
% idpar:     Vector of initial condition parameters
% vtype:     Selects potential type
% vpar:      Vector of potential parameters
% nLevels:   Number of levels for convergence test, must be 3 minimum
%
% Outputs
% t:         Vector of coarse t coordinates        [nt]
% scaledErr: Array of computed scaled dpsi norms   [nl x nt]

% determine number of levels used
num_lvls = lmax - lmin + 1;
% vector of range of levels
levels = lmin: lmax;

% solve lowest level sch_1d 
[x, t, ~, ~, ~, ~, ~, ~] = sch_1d_cn(tmax, lmin, lambda, idtype, idpar, vtype, vpar);

% solution vector for the scaled norms
scaledErr = zeros(num_lvls, length(t));

% getting parameter for exact solution for sch_1d
m = idpar(1);
% solve for exact solution 
exactSoln = exp(-1i * m^2 * pi^2 * t') .* sin(m * pi * x);

% iterate over all levels
for ls = 1: num_lvls
    % current level is
    level = levels(ls);

    % solve sch_1d at this level
    [~, ~, psi, ~, ~, ~, ~, ~] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);
    % match up correct mesh spacing
    psi = psi(1:2^(ls-1):end, 1:2^(ls-1):end);

    % calculate 2 norm of the error and scale error appropriately
    scaledErr(ls, :) = 4^(ls-1) * vecnorm(psi - exactSoln, 2, 2);
end
end