function [levels, t, scaledNorm] = sch_1d_lvl_err(tmax, lmin, lmax, lambda, idtype, idpar, vtype, vpar)
% Calculates the error between levels of the 1D Schrodinger Eq.
% 
% Inputs
% tmax:       max integration time
% lmin:       Min level, specifies the range over which to find errors
% lmax:       Max level, specifies the range over which to find errors
% lambda:     dt/dx
% idtype:     Selects initial condition type
% idpar:      Vector of initial condition parameters
% vtype:      Selects potential type
% vpar:       Vector of potential parameters
% nLevels:    Number of levels for convergence testing
%
% Outputs
% t:          Time coordinates for lowest level      [nt]
% scaledNorm: Array of scaled errors between levels  [nl x nt]

% determine number of levels to use
num_lvls = lmax - lmin;
% array of levels (subtract 1 because we calculate between levels)
levels = lmin: lmax-1;

% solve lowest level sch_1d and get time coordinates
[~, t, lowerLevelPsi, ~, ~, ~, ~, ~] = sch_1d_cn(tmax, lmin, lambda, idtype, idpar, vtype, vpar);
% vector to hold the norms. each row is a level
scaledNorm = zeros(num_lvls, length(t));

% iterate over the levels
for c_index = 1: num_lvls
    % current level, we skip over the lowest level (already solved)
    level = lmin + c_index;
    % find the solution to the 1d Schrodinger eq. for this level
    [~, ~, currentPsi, ~, ~, ~, ~, ~] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);
    % match up mesh spacing between levels in space/time
    currentPsi = currentPsi(1:2^c_index:end, 1:2^c_index:end);

    % compute the RMS of the difference between levels
    % choosing p=2 and dim=2 evaluates the norm row-wise (at each time)
    % apply scaling to errors (4^0, 4^1, 4^2, ...)
    scaledNorm(c_index, :) = 4^(c_index-1) * vecnorm(currentPsi - lowerLevelPsi, 2, 2);

    % update what the previous level is, and repeat
    lowerLevelPsi = currentPsi;
end
end