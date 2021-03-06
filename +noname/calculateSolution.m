% Calculates the solution of discretization for a given set of ms ts and orders.
%    discrHand -- function handle to a Discretization constructor
%    m         -- grid parameter
%    order     -- order of accuracy of the approximtion
%    T         -- time to calculate solution for
%    tsOpt     -- options for the time stepper creation.
%    input paramters m, t, order may all be vectors.
function [] = calculateSolution(filename, name, discrHand, m, T_in, order, tsOpt, force_flag)
    default_arg('force_flag',false);
    default_arg('tsOpt', []);

    if exist(filename,'file') && ~force_flag
        fprintf('File ''%s'' already exist.',filename);
        do_append = yesnoQuestion('Do you want to append to it?');
        if ~do_append
            fprintf('Exiting. No Solutions calculated.\n');
            return
        end
    end

    sf = SolutionFile(filename);

    orderWidth = findFieldWidth('%d',order);
    mWidth = findFieldWidth('%d',m);
    TWidth = findFieldWidth('%d',T_in);

    for i = 1:length(order)
        for j = 1:length(m)
            T = sort(T_in); % Make sure times are sorted

            discr = discrHand(m(j),order(i));
            k_max = discr.getTimestep(tsOpt);

            % Do we want to to save the initial conditions?
            if T(1) == 0
                snapshot = discr.getTimeSnapshot(0);
                saveToFile(sf, name, order(i), m(j),T(1), snapshot, NaN, NaN, discr);
                T(1) = [];
            end

            % Find out if times to be calulated are integer multiples of the smallest one.
            time_multiples = T/T(1);

            is_int_multiples = all(time_multiples == int64(time_multiples));

            if is_int_multiples
                fprintf('Calculating time series in increments\n');
            else
                fprintf('RESTARTING for each time in timeseries\n');
                fprintf('If this is not what you want try giving T in integer multiples.\n');
            end

            % T now contains all the times we need to step to,
            % if T contained 0 it has now been removed.

            if is_int_multiples
                % Times are integer multiples, we can save time
                [k,N] = alignedTimestep(k_max,T(1));
                tsOpt.k = k;
                ts = discr.getTimestepper(tsOpt);
                runtime = 0;
                for l = 1:length(T)
                    end_step = N * time_multiples(l);
                    fprintf('[order = %-*d, m = %-*d, T = %-*d]: ',orderWidth,order(i),mWidth,m(j),TWidth,T(l));
                    clock_start = tic();
                    ts.stepN(end_step-ts.n,true);
                    runtime = runtime + toc(clock_start);
                    snapshot = discr.getTimeSnapshot(ts);
                    saveToFile(sf, name, order(i), m(j),T(l), snapshot, runtime, k, discr);
                    fprintf('Done! (%.3fs)\n',runtime);
                end
            else
                % Times are not interger multiples, we have to start from 0 every time.
                for l = 1:length(T)
                    [k,N] = alignedTimestep(k_max,T(l));
                    tsOpt.k = k;
                    ts = discr.getTimestepper(tsOpt);
                    fprintf('[order = %-*d, m = %-*d, T = %-*d]: ',orderWidth,order(i),mWidth,m(j),TWidth,T(l));
                    clock_start = tic();
                    [v,t] = ts.stepN(N-ts.n,true);
                    runtime = toc(clock_start);
                    snapshot = discr.getTimeSnapshot(ts);
                    saveToFile(sf, name, order(i), m(j),T(l), snapshot, runtime, k, discr);
                    fprintf('Done! (%.3fs)\n',runtime);
                end

            end
            sf.stupidSave();
        end
    end
end


function saveToFile(sf, name, order, m, T, snapshot, runtime, k, discr)
    key.name  = name;
    key.order = order;
    key.m     = m;
    key.T     = T;

    entry.repr = snapshot;
    entry.runtime = runtime;
    entry.k = k;
    entry.discr = discr;

    sf.store(key,entry);
end
