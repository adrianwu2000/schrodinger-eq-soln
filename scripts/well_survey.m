% % well_survey.m determines how likely a particle is to be within a given
% range, with a well potential, and plots the excess fractional probability

% we want uniformly spaced points of ln(v)
% then we take exp(ln(v)) to get back v, the potential
% determines values for well depth
v0 = -exp(linspace(2, 10, 251));

% solution vector for excess fractional probability
F_e = zeros(size(v0));

% parameters
tmax = 0.10;
% level = 9;
level = 9;
lambda = 0.01;
idtype = 1; %(boosted Gaussian, p = 0)
idpar = [0.40, 0.075, 0];
vtype = 1;  %(rectangular well)
x1 = 0.6;
x2 = 0.8;

% iterate for each potential energy
for v = 1: length(v0)
    % solve 1D eq. for the probability and spatial mesh
    % depth of well is determined by values of v
    [x, ~, ~, ~, ~, ~, prob, ~] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, [0.6, 0.8, v0(v)]);
    nx = length(x);

    % find index where x1 and x2 are
    min = round(x1 * nx);
    max = round(x2 * nx);

    % compute temporal average
    temp_avg = mean(prob, 1);

    % normalize
    temp_avg = temp_avg / temp_avg(end);

    % compute excess fractional probability
    F_e(1, v) = (temp_avg(max) - temp_avg(min)) / (x2-x1);

end

% plotting
h = figure;
plot(log(v0), log(F_e));

title('Excess Probability vs. Well Potential')
subtitle('$$x_1 = 0.6,\quad x_2 = 0.8,\quad \ln(|V_0|) \in [2, 10]$$','interpreter','latex')
xlabel('$\ln(|V_0|)$', 'Interpreter', 'Latex');
ylabel('$\ln(\bar{F}_e(x_1, x_2))$', 'Interpreter', 'Latex');

% save figure
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'well_survey','-dpdf','-r0')
