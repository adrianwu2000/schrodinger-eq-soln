% ctest_1d.m performs convergence testing of sch_1d_cn.m for different
% cases as well as plots the results.

%%%%%%%%%%%%%%%%%%%%%%% EXACT FAMILY: EXACT ERROR %%%%%%%%%%%%%%%%%%%%%%%
clf;
% parameters: 
idtype = 0;
idpar = 3;
vtype = 0;
vpar = [];
tmax = 0.25; 
lambda = 0.1;
lmin = 6;
lmax = 9;

% calculate the error between the exact solution and each level
[levels, ts, scaledErrs] = sch_1d_exact_err(tmax, lmin, lmax, lambda, idtype, idpar, vtype, vpar);

h1 = figure(1);
% plot each level with the corresponding label
for ls = 1: length(levels)
    hold on;
    thisLabel = ['4^', num2str(ls-1), '||E(\psi^{', num2str(levels(ls)), '})||'];
    plot(ts, scaledErrs(ls, :), 'DisplayName', thisLabel);
end

%
xlabel('$t$', 'Interpreter', 'Latex');
ylabel('Scaled Errors');
title('Exact Family - Exact error');
legend show;

% save figure
set(h1,'Units','Inches');
pos = get(h1,'Position');
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h1,'exact_family_exact_err','-dpdf','-r0')

%%%%%%%%%%%%%%%%%%%%%%% EXACT FAMILY: LEVEL ERROR %%%%%%%%%%%%%%%%%%%%%%%

h2 = figure(2);
% calcualte the error between levels
[levels, tConv, lvlErr] = sch_1d_lvl_err(tmax, lmin, lmax, lambda, idtype, idpar, vtype, vpar);

% plot each level with the corresponding label
for ls = 1: length(levels)
    hold on;
    thisLabel = ['4^', num2str(ls-1), '||d\psi^{', num2str(levels(ls)), '}||'];
    plot(tConv, lvlErr(ls, :), 'DisplayName', thisLabel);
end

xlabel('$t$', 'Interpreter', 'Latex');
ylabel('Scaled Errors');
title('Exact Family - Error between levels');
legend show;

% save figure
set(h2,'Units','Inches');
pos = get(h2,'Position');
set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h2,'exact_family_lvl_err','-dpdf','-r0')

%%%%%%%%%%%%%%%%%%%%% BOOSTED GAUSSIAN: LEVEL ERROR %%%%%%%%%%%%%%%%%%%%%%
% parameters: 
idtype = 1;
idpar = [0.50 0.075 0.0];
vtype = 0;
vpar = [];
tmax = 0.01; 
lambda = 0.01;
lmin = 6;
lmax = 9;

h3 = figure(3);

% calculate the error between levels
[levels, tConv, boostLvlErr] = sch_1d_lvl_err(tmax, lmin, lmax, lambda, idtype, idpar, vtype, vpar);

% plot each level with the corresponding label
for ls = 1: length(levels)
    hold on;
    thisLabel = ['4^', num2str(ls-1), '||d\psi^{', num2str(levels(ls)), '}||'];
    plot(tConv, boostLvlErr(ls, :), 'DisplayName', thisLabel);
end

xlabel('$t$', 'Interpreter', 'Latex');
ylabel('Scaled Errors');
title('Boosted Gaussian - Error between levels');
legend show;

% save figure
set(h3,'Units','Inches');
pos = get(h3,'Position');
set(h3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h3,'boost_gauss_lvl_err','-dpdf','-r0')

