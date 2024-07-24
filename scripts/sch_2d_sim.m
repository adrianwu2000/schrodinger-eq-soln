% simulation for 2D Schrodinger Eq. over some situations (double/single
% slit, barrier/well) and make an animation
clearvars;

% choose experiment:
experiment = 5;
% 1: rectangular barrier (lower energy)
% 2: rectangular barrier (higher energy)
% 3: rectangular well    (shallower well)
% 4: rectangular well    (deeper well)
% 5: single slit
% 6: double slit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make plot
titles = {'Rectangular Barrier', 'Rectangular Barrier', 'Rectangular Well', 'Rectangular Well', 'Single Slit', 'Double Slit'};
fig = figure('Name', titles{experiment});
hold on;

% animate the chosen experiment
switch experiment
    case 1
        % rectangular barrier (low energy barrier)
        % parameters
        tmax = 0.04;
        level = 7;
        lambda = 0.01;
        idtype = 1;
        idpar = [0.25, 0.5, 0.1, 0.1, 10, 0];
        vtype = 1;
        vpar = [0.45, 0.55, 0, 1, 1E3];

        % animate using these parameters
        name = 'sim_2d_low_barrier.avi';
        sch_2d_animate(tmax, level, lambda, idtype, idpar, vtype, vpar, name);

    case 2
        % rectangular barrier (higher energy barrier)
        % parameters
        tmax = 0.04;
        level = 7;
        lambda = 0.01;
        idtype = 1;
        idpar = [0.25, 0.5, 0.1, 0.1, 10, 0];
        vtype = 1;
        vpar = [0.45, 0.55, 0, 1, 1E6];

        % animate using these parameters
        name = 'sim_2d_high_barrier.avi';
        sch_2d_animate(tmax, level, lambda, idtype, idpar, vtype, vpar, name);
    
    case 3
        % rectangular well (shallower energy well)
        % parameters
        tmax = 0.04;
        level = 7;
        lambda = 0.01;
        idtype = 1;
        idpar = [0.25, 0.5, 0.1, 0.1, 10, 0];
        vtype = 1;
        vpar = [0.45, 0.65, 0, 1, -1E3];

        % animate using these parameters
        name = 'sim_2d_shallow_well.avi';
        sch_2d_animate(tmax, level, lambda, idtype, idpar, vtype, vpar, name);

    case 4
        % rectangular well (deeper energy well)
        % parameters
        tmax = 0.04;
        level = 7;
        lambda = 0.01;
        idtype = 1;
        idpar = [0.25, 0.5, 0.1, 0.1, 10, 0];
        vtype = 1;
        vpar = [0.45, 0.65, 0, 0.8, -1E5];

        % animate using these parameters
        name = 'sim_2d_deep_well.avi';
        sch_2d_animate(tmax, level, lambda, idtype, idpar, vtype, vpar, name);

    case 5
        % single slit 
        % parameters
        tmax = 0.007;
        level = 6;
        lambda = 0.001;
        idtype = 1;
        idpar = [0.5, 0.125, 0.3, 0.05, 0, 100];
        vtype = 2;
        vpar = [0.45, 0.55, 0, 0, 1E6];

        % animate using these parameters
        name = 'sim_2d_single_slit.avi';
        sch_2d_animate(tmax, level, lambda, idtype, idpar, vtype, vpar, name);

    case 6
        % double slit
        % parameters
        tmax = 0.005;
        level = 7;
        lambda = 0.001;
        idtype = 1;
        idpar = [0.5, 0.125, 0.3, 0.05, 0, 100];
        vtype = 2;
        vpar = [0.35, 0.4, 0.6, 0.65, 1E6];

        % animate using these parameters
        name = 'sim_2d_double_slit.avi';
        sch_2d_animate(tmax, level, lambda, idtype, idpar, vtype, vpar, name);

end
