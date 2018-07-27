function A = spdiagsVariablePeriodic(vals,diags)
    % Creates an m x m periodic discretization matrix.
    % vals - m x ndiags matrix of values
    % diags - 1 x ndiags vector of the 'center diagonals' that vals end up on
    % vals that are not on main diagonal are going to spill over to
    % off-diagonal corners.

    default_arg('diags',0);

    [m, ~] = size(vals);

    A = sparse(m,m);

    for i = 1:length(diags)

        d = diags(i);
        a = vals(:,i);

        % Sub-diagonals
        if d < 0
            a_bulk = a(1+abs(d):end);
            a_corner = a(1:1+abs(d)-1);
            corner_diag = m-abs(d);
            A = A + spdiagVariable(a_bulk, d);
            A = A + spdiagVariable(a_corner, corner_diag);

        % Super-diagonals
        elseif d > 0
            a_bulk = a(1:end-d);
            a_corner = a(end-d+1:end);
            corner_diag = -m + d;
            A = A + spdiagVariable(a_bulk, d);
            A = A + spdiagVariable(a_corner, corner_diag);

        % Main diagonal
        else
             A = A + spdiagVariable(a, 0);
        end

    end

end