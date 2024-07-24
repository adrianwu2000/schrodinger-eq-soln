function sch_2d_animate(tmax, level, lambda, idtype, idpar, vtype, vpar, video_name)
% Animates the evolution of a solution to the 2D Schrodinger equation
% 
% Inputs
% tmax:       Maximum integration time
% level:      Discretization level
% lambda:     dt/dx
% idtype:     Selects initial condition type
% idpar:      Vector of initial condition parameters
% vtype:      Selects potential type
% vpar:       Vector of potential parameters
% video_name: Filename for the video to be saved as

    % get solution
    [x, y, t, ~, ~, ~, psimod, ~] = sch_2d_adi(tmax, level,lambda, idtype, idpar, vtype, vpar);
    num_steps = length(t);

    % plot and animate
    writerObj = VideoWriter(video_name);
    open(writerObj);
    
    % get potential type and get the subtitle
    switch vtype
        case 0
            % zero potential
            subt = 'Zero potential';
        case 1
            % well/barrier
            xmin = string(vpar(1));
            xmax = string(vpar(2));
            ymin = string(vpar(3));
            ymax = string(vpar(4));
            vc = string(vpar(5));
            subt = append('$$V_c = ', vc, ',\ \mathrm{for}\ x \in [', xmin, ', ', xmax, '], \ y \in [', ymin, ', ', ymax, ']$$');

        case 2
            % double/single slit
            x1 = string(vpar(1));
            x2 = string(vpar(2));
            x3 = string(vpar(3));
            x4 = string(vpar(4));
            vc = string(vpar(5));
            subt = append('$$\mathrm{Slits}\ \mathrm{located}\ \mathrm{at}\ [', x1, '\leq x \leq ', x2, ']\ \mathrm{and}\ [', x3, '\leq x \leq', x4, ']$$');
    end

    xlabel('$x$', 'Interpreter', 'Latex');
    ylabel('$y$', 'Interpreter', 'Latex');
    subtitle(subt,'interpreter','latex')
    % color scheme
    c = colorbar;
    c.Label.Interpreter = 'Latex';
    c.Label.String = '$|\Psi|$';
    clim([0, 1]);


    % animate each time step
    for t_step = 1: num_steps
        cla;

        % make contour plot
        contourf(x, y, squeeze(psimod(t_step, :, :)));

        % make color map have defined sections for colors
        C=parula(8);
        C(end+1,:)=1;
        colormap(flipud(C))

        % plot commands
        title(sprintf('Step %d of %d', t_step, num_steps));
        drawnow;
        writeVideo(writerObj, getframe(gcf));
        pause(0);
    end
    close(writerObj);
end