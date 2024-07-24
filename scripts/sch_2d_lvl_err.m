function [levels, t, scaledNorm] = sch_2d_lvl_err(tmax, lmin, lmax, lambda, idtype, idpar, vtype, vpar)
% Calculates the error between levels of the 1D Schrodinger Eq.
% 
% Inputs
% tmax:       max integration time
% lmin:       min level
% lmax:       max level
% lambda:     dt/dx
% idtype:     Selects initial condition type
% idpar:      Vector of initial condition parameters
% vtype:      Selects potential type
% vpar:       Vector of potential parameters
% nLevels:    Number of levels for convergence test (assumes at least 3)
%
% Outputs
% t:          Time coordinates for lowest level     [nt]
% scaledNorm: Array of scaled errors between levels [nl x nt]

% determine number of levels to use
num_lvls = lmax - lmin;
% the range of levels
levels = lmin: lmax-1;

% calculate solution of lowest level
[~, ~, t, lowerLevelPsi, ~, ~, ~, ~] = sch_2d_adi(tmax, lmin, lambda, idtype, idpar, vtype, vpar);

% solution vector for scaled errors
scaledNorm = zeros(num_lvls, length(t));

for ls = 1: num_lvls
    % the current level is
    level = lmin + ls;
    % solution at current level
    [~, ~, ~, currentPsi, ~, ~, ~, ~] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);
    % match up mesh spacing correctly based on level
    currentPsi = currentPsi(1:2^ls:end, 1:2^ls:end, 1:2^ls:end);
    % calculate root mean square of the error between levels
    % we use vecnorm() twice because there are two directions (nx, ny) to
    % sum over
    dpsiNorm = vecnorm(vecnorm(currentPsi - lowerLevelPsi, 2, 2), 2, 3);
    % scale appropriately (4^0, 4^1, 4^2, ...)
    scaledNorm(ls, :) = 4^(ls-1) * dpsiNorm;
    % update the lower level and repeat process
    lowerLevelPsi = currentPsi;
end
end