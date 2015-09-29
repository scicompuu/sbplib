% Reference is either a key or a function handle
function [q, e, h] = convergence(filename, errorFunc, reference, method, order, m, T)
    default_arg('errorFunc', @scheme.error1d);

    sf = SolutionFile(filename);

    analytical_ref = isa(reference,'function_handle');
    if ~analytical_ref
        reference = sf.get(reference);
    end


    % Generate convergence, error, and efficiency plots for each search key with more than one entry.
    for i = 1:length(m)
        key.method = method;
        key.order = order;
        key.m = m(i);
        key.T = T;

        entry = sf.get(key);

        [e(i),h(i)] = errorForEntry(key, entry, errorFunc, reference,T);

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
        x = v_repr.x;
        v_ref = reference(x,T);
    else
        % Downsample the reference solution
        x = v_repr.x;
        x_ref = reference.x;

        [~,I] = ismember(x,x_ref,'rows');
        v_ref = reference.v(I);
    end

    e = errorFunc(discr,v, v_ref);
    h = discr.h;
end