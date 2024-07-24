function [x, y, t, psi, psire, psiim, psimod, v] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar)
% Inputs
%
% tmax:   Maximum integration time
% level:  Discretization level
% lambda: dt/dx
% idtype: Selects initial data type
% idpar:  Vector of initial data parameters
% vtype:  Selects potential type
% vpar:   Vector of potential parameters
%
% Outputs
%
% x:      Vector of x coordinates                  [nx]
% y:      Vector of y coordinates                  [ny]
% t:      Vector of t coordinates                  [nt]
% psi:    Array of computed psi values             [nt x nx x ny]
% psire:  Array of computed psi_re values          [nt x nx x ny]
% psiim:  Array of computed psi_im values          [nt x nx x ny]
% psimod: Array of computed sqrt(psi psi*) values  [nt x nx x ny]
% v:      Array of potential values                [nx x ny]

% get discretized time/space domain while keeping ratio constant
% copied from assignment description
nx = 2^level + 1;
x = linspace(0.0, 1.0, nx);
y = linspace(0.0, 1.0, nx);      
dx = x(2) - x(1); % dx = dy

% create time coordinates 
dt = lambda * dx;
nt = round(tmax / dt) + 1;
t = [0:nt-1] * dt;

% initialize psi
psi = zeros(nt, nx, nx);

% determine initial condition
if idtype == 0
    % parameters
    m_x = idpar(1);
    m_y = idpar(2);
    % exact family
    psi(1, :, :) = sin(m_x * pi * x) .* sin(m_y * pi * y');

elseif idtype == 1
    % parameters
    x_0 = idpar(1);
    y_0 = idpar(2);
    delta_x = idpar(3);
    delta_y = idpar(4);
    p_x = idpar(5);
    p_y = idpar(6);
    % boosted gaussian
    psi(1, :, :) = exp(1i * (p_x * x + p_y * y') - ((x - x_0) ./ delta_x) .^ 2 - ((y' - y_0) ./ delta_y) .^ 2);
end

% determine potential

% case when vtype == 0, zero potential
v = zeros(nx, nx);

if vtype == 1
    % parameters 
    xmin = vpar(1);
    xmax = vpar(2);
    ymin = vpar(3);
    ymax = vpar(4);
    vc = vpar(5);

    % rectangular barrier/well
    % find index where min/max are
    ymin_index = round(ymin * (nx-1)) + 1;
    ymax_index = round(ymax * (nx-1)) + 1;
    xmin_index = round(xmin * (nx-1)) + 1;
    xmax_index = round(xmax * (nx-1)) + 1;

    % change potential to correct value in that range
    v(ymin_index:ymax_index, xmin_index:xmax_index) = vc;

elseif vtype == 2
    % parameters
    x1 = vpar(1);
    x2 = vpar(2);
    x3 = vpar(3);
    x4 = vpar(4);
    vc = vpar(5);

    % double slit
    % find index where slits are
    x1_i = round(x1 * (nx-1)) + 1;
    x2_i = round(x2 * (nx-1)) + 1;
    x3_i = round(x3 * (nx-1)) + 1;
    x4_i = round(x4 * (nx-1)) + 1;

    % change potential to correct value
    j_prime = 1 + (nx - 1) / 4;

    % wall at location j_prime (thickness 2)
    v(j_prime: j_prime+1, :) = vc;

    % slits along that wall between (x1,x2) and (x3,x4)
    v(j_prime: j_prime+1, x1_i:x2_i) = 0;
    v(j_prime: j_prime+1, x3_i:x4_i) = 0;
end

% set up tridiagonal system

% Initialize storage for sparse matrix and RHS
dl = zeros(nx,1);
d  = zeros(nx,1);
du = zeros(nx,1);
f  = zeros(nx,1); % RHS
% B = zeros(nx,1);

% define LHS coefficients
a = dt / dx^2;
dl = -0.5i * a * ones(nx, 1);
d = (1 + 1i * a) * ones(nx, 1);
du = dl;
% boundary cases
d([1, nx]) = 1.0; % edge
dl(nx-1) = 0.0;   % inner
du(2) = 0.0;      % inner

% make tridiagonal system
A = spdiags([dl, d, du], -1:1, nx, nx);

% use ADI to solve 2D Schrodinger eq.
% iterate over all time steps
for tstep = 1: nt-1
    % step 1
    for jj = 2:nx-1
        % For each (j), we solve a tridiagonal system of the intermediate
        % grid function for all columns of (i).

        % Define RHS of the linear system
        % --------------------------------------------------|---|---|---|-
        p_ij = psi(tstep, 2:nx-1, jj); %                    | n | i | j |
        % --------------------------------------------------|---|---|---|-
        p_ip1j = psi(tstep, 3:nx, jj); %                    | n |i+1| j |
        % --------------------------------------------------|---|---|---|-
        p_im1j = psi(tstep, 1:nx-2, jj); %                  | n |i-1| j | 
        % --------------------------------------------------|---|---|---|-
        p_ijm1 = psi(tstep, 2:nx-1, jj-1); %                | n | i |j-1|
        % --------------------------------------------------|---|---|---|-
        p_ijp1 = psi(tstep, 2:nx-1, jj+1); %                | n | i |j+1|
        % --------------------------------------------------|---|---|---|-
        p_im1jm1 = psi(tstep, 1:nx-2, jj-1); %              | n |i-1|j-1|
        % --------------------------------------------------|---|---|---|-
        p_ip1jm1 = psi(tstep, 3:nx, jj-1); %                | n |i+1|j-1|
        % --------------------------------------------------|---|---|---|-
        p_im1jp1 = psi(tstep, 1:nx-2, jj+1); %              | n |i-1|j+1|
        % --------------------------------------------------|---|---|---|-
        p_ip1jp1 = psi(tstep, 3:nx, jj+1); %                | n |i+1|j+1|
        % --------------------------------------------------|---|---|---|-
        v_ij = v(2:nx-1, jj)'; %                            | n | i | j |
        % --------------------------------------------------|---|---|---|-

        % constants
        c1 = (1 - 1i * a) * ((1 - 1i * a) - 0.5i * dt * v_ij);
        c2 = 0.5i * a * (1 - 1i * a - 0.5i * dt * v_ij);
        c3 = 0.5i * a * (1 - 1i * a);
        c4 = -0.25 * a * a;
        
        % plug and solve
        f(2:nx-1) = c1 .* p_ij + c2 .* (p_im1j + p_ip1j) + c3 * (p_ijm1 + p_ijp1) + c4 * (p_im1jm1 + p_ip1jm1 + p_im1jp1 + p_ip1jp1);
        % 
        % %%%%%%%%%%%%%%%%%%% test %%%%%%%%%%%%%%%%%%%%%%%%%
        % bracket = (p_ijp1 - 2 *p_ij + p_ijm1) / (dx*2);
        % B(2:nx-1) = p_ij + 1i * dt/2 * (bracket) - 1i * (dt/2) * v_ij' * p_ij;
        % 
        % B_ij = B(tstep, 2:nx-1, jj);     % A^n_i, j
        % B_ip1j = B(tstep, 3:nx, jj);     % A^n_i+1, j
        % B_im1j= B(tstep, 1:nx-2, jj);    % A^n_i-1, j
        % 
        % f_test(2:nx-1) = B_ij + 1i * (dt/2) * ( (B_ip1j - 2 * B_ij + B_im1j) / (dx^2) );
        % 
        % 
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % solve for intermediete linear system
        psi(tstep+1, :, jj) = A \ f;
    end
    % step 2
    for ii = 2:nx-1
        % we flip - for each (i) we solve a tridiagonal system of the
        % intermediete grid function for all rows of (j)
        
        % edit tridiagonal system
        % du, dl are the same
        d = (1 + 1i * a) * ones(nx, 1) + (1/2)*1i * dt * v(ii, :)';
        d([1, nx]) = 1.0;
        % RHS
        B = spdiags([dl, d, du], -1:1, nx, nx);
        % solve system using intermediete solution
        psi(tstep+1, ii, :) = B \ squeeze(psi(tstep+1, ii, :));
    end
end
% outputs
psiim = imag(psi);
psire = real(psi);
psimod = abs(psi);
end