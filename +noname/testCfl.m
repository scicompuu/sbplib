% noname.testCfl(discr, timestepper_method, T, alpha0, tol,threshold, silentFlag)
% Example:
% noname.testCfl(Discr(100,4), 'rk4', 1, [0, 1])
function testCfl(discr, timestepper_method, T, alpha0, tol,threshold, silentFlag)
    default_arg('tol',0.00005);
    default_arg('threshold',1e2);
    default_arg('silentFlag', false);

    % TODO:
    % Set threshold from the initial conditions of the pde?
    % Take a set number of steps instead of evolving to a certain time?
    % Stop evolving when it has blown up?

    testAlpha = getAlphaTester(discr, T, threshold, silentFlag, timestepper_method);

    % Make sure that the upper bound is not working
    ok = testAlpha(alpha0(2));
    if ok % Upper bound too large!
        error('The upper bound on alpha is stable!')
    end

    % Make sure that the lower bound is ok
    if alpha0(1) ~= 0
        ok = testAlpha(alpha0(1));
        if ~ok
            error('The lower bound on alpha is unstable!');
        end
    end

    if silentFlag
        rsInterval = util.ReplaceableString('');
    end

    % Use bisection to find sharp estimate
    while( (alpha0(2)-alpha0(1))/alpha0(1) > tol)
        alpha = mean(alpha0);

        if ~silentFlag
            fprintf('[%.3e,%.3e]: ', alpha0(1), alpha0(2));
        else
            rsInterval.update('[%.3e,%.3e]: ', alpha0(1), alpha0(2));
        end

        [ok, n_step, maxVal] = testAlpha(alpha);

        if ok
            alpha0(1) = alpha;
            stability = 'STABLE';
        else
            alpha0(2) = alpha;
            stability = 'UNSTABLE';
        end

        if ~silentFlag
            fprintf('a = %.3e, n_step=%d %8s max = %.2e\n', alpha, n_step, stability, maxVal);
        end
    end

    if silentFlag
        rsInterval = util.ReplaceableString('');
    end

    fprintf('T = %-3d dof = %-4d order = %d: clf = %.4e\n',T, discr.size(), discr.order, alpha0(1));

end

function f = getAlphaTester(discr, T, threshold, silentFlag, timestepper_method)

    % Returns true if cfl was ok
    function [ok, n_step, maxVal] = testAlpha(alpha)
        ts = discr.getTimestepper(struct('method', timestepper_method, 'cfl', alpha));

        warning('off','all')
        ts.evolve(T,true);
        warning('on','all')

        [v,t] = ts.getV();
        maxVal = max(v);

        if isnan(maxVal) || maxVal == Inf || abs(maxVal) > threshold
            ok = false;
        else
            ok = true;
        end

        n_step = ts.n;
    end

    f = @testAlpha;
end
