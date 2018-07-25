% [discr, trueSolution] =  schemeFactory(m)
%     where trueSolution should be a timeSnapshot of the true solution a time T
% T is the end time
% m are grid size parameters.
% N are number of timesteps to use for each gird size
% timeOpt are options for the timeStepper
% errorFun is a function_handle taking 2 or 3 arguments, errorFun(trueSolution, approxSolution), errorFun(trueSolution, approxSolution, discr)
function e = calculateErrors(schemeFactory, T, m, N, errorFun, timeOpt)
    %TODO: Ability to choose paralell or not
    assertType(schemeFactory, 'function_handle');
    assertNumberOfArguments(schemeFactory, 1);
    assertScalar(T);
    assert(length(m) == length(N), 'Vectors m and N must have the same length');
    assertType(errorFun, 'function_handle');

    if ~ismember(nargin(errorFun), [2,3])
        error('sbplib:noname:calculateErrors:wrongNumberOfArguments', '"%s" must have 2 or 3, found %d', toString(errorFun), nargin(errorFun));
    end

    default_arg('timeOpt', struct());


    e = zeros(1,length(m));
    parfor i = 1:length(m)
        done = timeTask('m = %3d ', m(i));

        [discr, trueSolution] = schemeFactory(m(i));

        timeOptTemp = timeOpt;
        timeOptTemp.k = T/N(i);
        ts = discr.getTimestepper(timeOptTemp);
        ts.stepTo(N(i), true);
        approxSolution = discr.getTimeSnapshot(ts);

        switch nargin(errorFun)
            case 2
                e(i) = errorFun(trueSolution, approxSolution);
            case 3
                e(i) = errorFun(trueSolution, approxSolution, discr);
        end

        fprintf('e = %.4e', e(i))
        done()
    end
    fprintf('\n')
end


%% Example error function
% u_true = grid.evalOn(dr.grid, @(x,y)trueSolution(T,x,y));
% err = u_true-u_false;
% e(i) = norm(err)/norm(u_true);
% % e(i) = sqrt(err'*d.H*d.J*err/(u_true'*d.H*d.J*u_true));
