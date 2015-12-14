% noname.testCfl(discr, timestepper_method, T, alpha0, tol,threshold, silentFlag)
% Example:
% noname.testCfl(Discr(100,4), 'rk4', 1, [0, 1])
function testCfl(discr, timestepper_method, T, alpha0, tol,threshold, silentFlag)
    default_arg('tol',0.00005);
    default_arg('threshold',1e2);
    default_arg('silentFlag', false);

    if T < alpha0(2)
        error('Upper bound on alpha must be smaller than T');
    end

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

    % Use bisection to find sharp estimate
    while( (alpha0(2)-alpha0(1))/alpha0(1) > tol)
        alpha = mean(alpha0);
        fprintf('[%.3e,%.3e]: ', alpha0(1), alpha0(2));
        ok = testAlpha(alpha);
        if ok
            alpha0(1) = alpha;
        else
            alpha0(2) = alpha;
        end
    end

    fprintf('T = %-3d dof = %-4d order = %d: clf = %.4e\n',T, discr.size(), discr.order, alpha0(1));

end

function f = getAlphaTester(discr, T, threshold, silentFlag, timestepper_method)

    % Returns true if cfl was ok
    function ok = testAlpha(alpha)
        ts = discr.getTimestepper(struct('method', timestepper_method, 'cfl', alpha));

        warning('off','all')
        ts.evolve(T,true);
        warning('on','all')

        [v,t] = ts.getV();
        maxVal = max(v);

        if isnan(maxVal) || maxVal == Inf || maxVal > threshold
            ok = false;
        else
            ok = true;
        end

        if ~silentFlag
            fprintf('a = %.3e, max= %.2e\n', alpha, maxVal);
        end
    end

    f = @testAlpha;
end