% noname.animate(discretization, time_modifier, Tend, dirname, opt)
%
% Example:
%      noname.animate(discr,timemodifier,tend)
%      noname.animate(discr,1, [tstart tend],'my_mov', opt)

function animate(discretization, time_modifier, Tend, dirname, opt)
    default_arg('time_modifier', 1);
    default_arg('Tend', Inf);
    default_arg('dirname', '');

    optDefault.plotType = 'animation';
    optDefault.time = [];

    default_struct('opt', optDefault);


    if time_modifier < 0
        do_pause = true;
        time_modifier = -time_modifier;
    else
        do_pause = false;
    end

    if isinf(time_modifier)
        do_step = true;
    else
        do_step = false;
    end

    makemovies = ~strcmp(dirname,'');
    if makemovies
        dirname = ['mov/' dirname];
    end

    fprintf('Animating: %s\n',discretization.name);
    fprintf('order    : %d\n',discretization.order);
    fprintf('m        : %d\n',size(discretization));


    ts = discretization.getTimestepper(opt.time);

    if numel(Tend) == 2
        Tstart = Tend(1);
        Tend = Tend(2);


        fprintf('Evolving to starting time: ');
        ts.evolve(Tstart,'true');
        fprintf(' - Done\n');
        start_solution = discretization.getTimeSnapshot(ts);
    else
        start_solution = discretization.getTimeSnapshot(0);
        Tstart = start_solution.t;
    end

    [update, figure_handle] = discretization.setupPlot(opt.plotType);
    if makemovies
        save_frame = anim.setup_fig_mov(figure_handle,dirname);
    end

    % Initialize loop
    str = '';
    % Loop function
    function next_t = G(next_t)
        ts.evolve(next_t);
        sol = discretization.getTimeSnapshot(ts);
        update(sol);
        drawnow
        % waitforbuttonpress
        if makemovies
            save_frame();
        end

        if do_pause
            pause
        end
    end
    update(start_solution);

    fprintf('Using time step k = %.6f\n',ts.k);
    fprintf('System size: %d\n',size(discretization));
    % waitforbuttonpress


    if ~do_step
        pause
        anim.animate(@G, Tstart, Tend, time_modifier);
    else
        pause
        while ts.t < Tend
            ts.step();
            sol = discretization.getTimeSnapshot(ts);
            update(sol);
            drawnow

            if do_pause
                pause
            end
        end
    end

    % str = util.replace_string(str,'');

    % if makemovies
        % fprintf('Generating movies...\n')
        % system(sprintf('bash make_movie.sh %s',dirname));
    % end
end
