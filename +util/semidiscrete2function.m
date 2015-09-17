function F = semidiscrete2function(D,S)
    % TODO: Find threshhold when sparse makes sense.
    % TODO: Why does having the spase block slow down convergence.m (or does it?)
    if size(D,1) > 500
        D = sparse(D);
    end

    if ~isa(S,'function_handle')
        F = @(w,t)(D*w + S);
    else
        F = @(w,t)(D*w + S(t));
    end
end
