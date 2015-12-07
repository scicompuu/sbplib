function testCfl(discr, timestepper_method, T, alpha0, tol,threshold, silentFlag)
    default_arg('tol',0.00005);
    default_arg('threshold',1e2);
    default_arg('silentFlag', false);

    alpha0(2) = alpha0(1)+2*(alpha0(2)-alpha0(1));

    while( (alpha0(2)-alpha0(1))/alpha0(1) > tol)
        alpha = mean(alpha0);

        ts = discr.getTimestepper(struct('method', timestepper_method, 'cfl', alpha));

        warning('off','all')
        ts.evolve(T,true);
        warning('on','all')

        [v,t] = ts.getV();
        max_val = max(v);

        if isnan(max_val) || max_val == Inf || max_val > threshold
            alpha0(2) = alpha;
        else
            alpha0(1) = alpha;
        end

        if ~silentFlag
            fprintf('[%.3e,%.3e]: a = %.3e, max= %.2e\n',alpha0(1),alpha0(2), alpha,max_val);
        end
    end

    fprintf('T = %-3d dof = %-4d order = %d: clf = %.4e\n',T, discr.size(), discr.order, alpha0(1));

end