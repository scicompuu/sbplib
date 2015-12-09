% animate(dirname,discretization,Tend, time_modifier,time_method)
%
% Example:
%      animate('',discr,tend)
%      animate('my_mov',discr,tend,time_mod,time_method)

function hand = animate(discretization, time_modifier, Tend, dirname, opt)
    default_arg('time_modifier',1);
    default_arg('Tend', Inf);
    default_arg('dirname','');
    default_arg('opt', []);


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


    ts = discretization.getTimestepper(opt);

    if numel(Tend) == 2
        Tstart = Tend(1);
        Tend = Tend(2);


        fprintf('Evolving to starting time: ');
        ts.evolve(Tstart,'true');
        fprintf(' - Done\n');
        start_solution = discretization.getTimeSnapshot(ts);
    else
        Tstart = 0;
        start_solution = discretization.getTimeSnapshot(0);
    end

    [update, figure_handle] = discretization.setupPlot('animation');
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
        % pause(0.1)
        str = util.replace_string(str,'t = %.5f',ts.t);

        if do_pause
            pause
        end
    end
    update(start_solution);

    fprintf('Using time step k = %.6f\n',ts.k);
    fprintf('System size: %d\n',size(discretization));
    % waitforbuttonpress


    if ~do_step
        anim.animate(@G, Tstart, Tend, time_modifier);
    else
        while true
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
