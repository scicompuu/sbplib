% animate(dirname,discretization,Tend, time_modifier,time_method)
%
% Example:
%      animate('',discr,tend)
%      animate('my_mov',discr,tend,time_mod,time_method)

function hand = animate(dirname,discretization,Tend, time_modifier,time_method)
    makemovies = ~strcmp(dirname,'');
    if makemovies
        dirname = ['mov/' dirname];
    end

    default_arg('Tend',Inf);
    default_arg('time_modifier',1);
    default_arg('time_method',[]);

    if time_modifier < 0
        do_pause = true;
        time_modifier = -time_modifier;
    else
        do_pause = false;
    end


    fprintf('Animating: %s\n',discretization.name);
    fprintf('Tend     : %.2f\n',Tend);
    fprintf('order    : %d\n',discretization.order);
    fprintf('m        : %d\n',size(discretization));

    fprintf('Creating time discretization');
    tic
    ts = discretization.getTimestepper(time_method);
    fprintf(' - done  %fs\n', toc())

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
    sol = discretization.getTimeSnapshot(0);
    update(sol);

    fprintf('Using time step k = %.6f\n',ts.k)
    fprintf('System size: %d\n',size(discretization))
    waitforbuttonpress
    anim.animate(@G,0,Tend,time_modifier)
    str = util.replace_string(str,'');

    % if makemovies
        % fprintf('Generating movies...\n')
        % system(sprintf('bash make_movie.sh %s',dirname));
    % end
end
