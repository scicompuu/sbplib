% Reshapes a matrix as if it was stored in row major order.
function B = reshapeRowMaj(A, m)
    D = length(m);

    if D == 1
        m = [m 1];
        D = 2;
    end

    % Reshape and reverse order of indecies
    B = permute(reshape(permute(A, ndims(A):-1:1), rot90(m,2)), D:-1:1);
end