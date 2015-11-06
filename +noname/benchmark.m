% animate(discretization, N, time_method)
%
% Example:
%      benchmark(discr,100)
%      benchmark(discr,1000,'rk4')

function hand = benchmark(discretization,N ,time_method,do_profile)
    default_arg('N',100);
    default_arg('time_method',[]);

    fprintf('Creating time discretization');
    tic
    ts = discretization.getTimestepper(time_method);
    fprintf(' - done  %fs\n', toc());

    if do_profile
        profile on
    end

    fprintf('Taking %d steps',N);
    tic;
    ts.stepN(N,true);
    fprintf(' - done  %fs\n', toc());

    if do_profile
        profile viewer
    end
end
