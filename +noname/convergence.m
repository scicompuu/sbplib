% Reference is either a key or a function handle
function [q, e, h, runtime] = convergence(filename, errorFunc, reference, name, order, m, T)
    default_arg('errorFunc', @scheme.error1d);

    sf = SolutionFile(filename);


    % Generate convergence, error, and efficiency plots for each search key with more than one entry.
    for i = 1:length(m)
        key.name = name;
        key.order = order;
        key.m = m(i);
        key.T = T;

        entry = sf.get(key);

        [e(i),h(i)] = errorForEntry(key, entry, errorFunc, reference,T);
        runtime(i) = entry.runtime;

    end
    q = convergence(e,h);
end

function [e, h] = errorForEntry(key,entry, errorFunc, reference,T)
    v_repr = entry.repr;
    discr = entry.discr;

    % Get the solution to be compared
    v = v_repr.v;

    % Get the reference solution vector
    if isa(reference,'function_handle');
        x = v_repr.grid.points();
        v_ref = reference(x,T);
    else
        % Downsample the reference solution
        v_ref = reference.grid.restrictFunc(reference.v, v_repr.grid);
    end

    e = errorFunc(discr,v, v_ref);
    h = discr.h;
end