% Calculates the solution of a discr at a given time using aligned timesteps.
% Returns the solution as a grid function as defined in +grid
function gf = calculateSolution(discr, T, tsOpt)
    k_max = discr.getTimestep(tsOpt);
    [k,N] = alignedTimestep(k_max,T);
    tsOpt.k = k;
    ts = discr.getTimestepper(tsOpt);

    gf = ts.stepN(N-ts.n);
end