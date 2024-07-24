% Convergenve test for 2D Schrodinger Eq. with Exact family initial
% conditions and zero potential
clf;

% parameters
idtype = 0;
vtype = 0;
idpar = [2, 3];
vpar = [];
tmax = 0.05;
lambda = 0.05;
lmin = 6;
lmax = 9;

%%%%%%%%%%%%%%%%%%%%%%% EXACT FAMILY: LEVEL ERROR %%%%%%%%%%%%%%%%%%%%%%%
h1 = figure(1);

% convergence test
[levels, t, scaledNorms] = sch_2d_lvl_err(tmax, lmin, lmax, lambda, idtype, idpar, vtype, vpar);

% plot each set of errors
for ls = 1: length(levels)
    hold on;
    label = ['4^', num2str(ls-1), '||d\psi^{', num2str(levels(ls)), '}||'];
    plot(t, scaledNorms(ls, :), 'DisplayName', label);
end
% plot
title('Exact Family - Error between levels');
subtitle('$$\ell = 6, 7, 8, 9$$','interpreter','latex')
xlabel('$t$', 'Interpreter', 'Latex');
ylabel('Scaled Errors')
legend show;

% save figure
set(h1,'Units','Inches');
pos = get(h1,'Position');
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h1,'2d_exact_fam_lvl_err','-dpdf','-r0')

%%%%%%%%%%%%%%%%%%%%%%% EXACT FAMILY: EXACT ERROR %%%%%%%%%%%%%%%%%%%%%%%
h2 = figure(2);

% convergence test
[levels, t, scaledNorm] = sch_2d_exact_err(tmax, lmin, lmax, lambda, idtype, idpar, vtype, vpar);

% plot each set of errors
for ls = 1: length(levels)
    hold on;
    label = ['4^', num2str(ls-1), '||d\psi^{', num2str(levels(ls)), '}||'];
    plot(t, scaledNorm(ls, :), 'DisplayName', label);
end

% plot
title('Exact Family - Exact error');
subtitle('$$\ell = 6, 7, 8, 9$$','interpreter','latex')
xlabel('$t$', 'Interpreter', 'Latex');
ylabel('Scaled Errors')
legend show;

% save figure
set(h2,'Units','Inches');
pos = get(h2,'Position');
set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h2,'2d_exact_family_exact_err','-dpdf','-r0')