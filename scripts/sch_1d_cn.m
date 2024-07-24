function [x, t, psi, psire, psiim, psimod, prob, v] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar)
% Solve the 1D SE problem using the Crank-Nicolson scheme.
% 
% Inputs
% tmax:   Maximum integration time
% level:  Discretization level
% lambda: dt/dx
% idtype: Selects initial xfrog_assets type
% idpar:  Vector of initial xfrog_assets parameters
% vtype:  Selects potential type
% vpar:   Vector of potential parameters
%
% Outputs
% x:      Vector of x coordinates                   [nx]
% t:      Vector of t coordinates                   [nt]
% psi:    Array of computed psi values              [nt x nx]
% psire:  Array of computed psi_re values           [nt x nx]
% psiim:  Array of computed psi_im values           [nt x nx]
% psimod: Array of computed sqrt(psi psi*) values   [nt x nx]
% prob:   Array of computed running integral values [nt x nx]
% v:      Array of potential values                 [nx]

% get discretized time/space domain while keeping ratio constant
% copied from diff_1d_imp.m
nx = 2^level + 1;
x = linspace(0.0, 1.0, nx);
dx = x(2) - x(1);
dt = lambda * dx;
nt = round(tmax / dt) + 1;
t = [0: nt-1] * dt;

% wave function psi solution vector for all times 
psi = zeros(nt, nx);

% check initial data type
% Set initial condition based on idtype
switch idtype
    case 0  % Exact family
        m = idpar(1);
        psi(1, :) = sin(m * pi * x);
    case 1  % Boosted Gaussian
        x0 = idpar(1);
        delta = idpar(2);
        rho = idpar(3);
        psi(1, :) = exp(1i * rho * x - ((x - x0) ./ delta) .^ 2);
end

% make potential V(x)
v = zeros(1, nx);
% check type of potential
switch vtype
    case 0
        % no change to V(x), still zero
    case 1
        % rectangular barrier / well parameters
        xmin = vpar(1);
        xmax = vpar(2);
        vc = vpar(3);
        % figure out what index xmin/xmax are
        min_index = round(xmin * (nx-1)) + 1;
        max_index = round(xmax * (nx-1)) + 1;
        % update
        v(min_index:max_index) = vc;
end

% set up tridiagonal system; copied from diff_1d_imp.m

% Initialize storage for sparse matrix and RHS
dl = zeros(nx,1); % c^+
d  = zeros(nx,1); % c^0
du = zeros(nx,1); % c^-
f  = zeros(nx,1); % RHS

% set up tridiagonal system
% LHS coefficients based on algebraic solution of Crank-Nicholson scheme
dl = 0.5 / dx^2 * ones(nx, 1);
d = (1i / dt - 1 / dx^2) * ones(nx, 1) - 0.5 * v';
du = dl;

% fix up boundary cases
d([1, nx]) = 1.0;   % first, last
du(2) = 0.0;        % inner point
dl(nx-1) = 0.0;     % inner point
% define sparse matrix
A = spdiags([dl, d, du], -1:1, nx, nx);

% compute solution using implicit scheme, iterate over all time steps
% disp('start loop')
% disp('total loops is: ')
% disp(nt-1)

for t_step = 1:nt-1

    % define RHS of linear system
    psi_n_j = psi(t_step, 2:nx-1);   % n, j
    psi_n_jm1 = psi(t_step, 1:nx-2); % n, j-1
    psi_n_jp1 = psi(t_step, 3:nx);   % n, j+1
    v1 = v(2:nx-1);                  % potential

    % plug and solve for RHS
    f(2:nx-1) = ((1i/dt) + 0.5 * v1) .* psi_n_j - (1/2) * (psi_n_jm1 - 2 * psi_n_j + psi_n_jp1) / dx^2;
    %f(1) = 0.0;
    %f(nx) = 0.0;

    % solve using left division
    psi(t_step+1, :) = A \ f;
end

% create outputs
psire = real(psi);
psiim = imag(psi);
psimod = abs(psi);
% psimod = sqrt(psi .* conj(psi));

% calcualte probability density running integral
prob = zeros(nt, nx);
% first step
prob(:, 1) = dx * psimod(:, 1) .^ 2;

for x_step = 2: nx
    % apply trapezoidal approximation for integral
    % multiply by 0.5?
    prob(:, x_step) = prob(:, x_step-1) + dx * psimod(:, x_step) .^ 2;
end

end