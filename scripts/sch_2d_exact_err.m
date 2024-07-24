function [levels, t, scaledNorm] = sch_2d_exact_err(tmax, lmin, lmax, lambda, idtype, idpar, vtype, vpar)
% Calculates the error between the 2D Schrodinger Eq. and the exact
% solution
% 
% Inputs
% tmax:       Maximum integration time
% level:      Discretization level
% lambda:     dt/dx
% idtype:     Selects initial condition type
% idpar:      Vector of initial condition parameters
% vtype:      Selects potential type
% vpar:       Vector of potential parameters
% nLevels:    Number of levels for convergence test, must be 3 minimum
%
% Outputs
% t:          Vector of coarse t coordinates        [nt]
% scaledNorm: Array of computed scaled dpsi norms   [nl x nt]

% determine number of levels to use
num_lvls = lmax - lmin + 1;
% the range of levels
levels = lmin: lmax;

% calculate coordinates for x,y,t of lowest level
[x, y, t, ~, ~, ~, ~, ~] = sch_2d_adi(tmax, lmin, lambda, idtype, idpar, vtype, vpar);
% solution array to store scaled errors
scaledNorm = zeros(num_lvls, length(t));
% assumes exact family initial condition
m_x = idpar(1);
m_y = idpar(2);
% calculate exact solution
psiTrue = exp(-1i * (m_x^2 + m_y^2) * pi^2 * t') .* sin(m_x * pi * permute(x, [1, 3, 2])) .* sin(m_y * pi * y);

% calculate error between exact solution and each level
for ls = 1: num_lvls
    % the current level is
    level = lmin + ls - 1;
    % calculate solution at this level
    [~, ~, ~, psi, ~, ~, ~, ~] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);
    % match up mesh grid spacing according to level
    psi = psi(1:2^(ls-1):end, 1:2^(ls-1):end, 1:2^(ls-1):end);
    % root mean square of error 
    unscaledNorm = vecnorm(vecnorm(psi - psiTrue, 2, 2), 2, 3);
    % scale error (4^0, 4^1, 4^2, ...)
    scaledNorm(ls, :) = 4^(ls-1) * unscaledNorm;
end
end