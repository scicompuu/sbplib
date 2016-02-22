% Takes a grid function and reshapes it into a matrix of shape m.
function F = reshapeKronVector(gf, m)
    D = length(m);

    if D == 1
        F = gf;
        return
    end

    % Reshape and reverse order of indecies
    F = permute(reshape(gf, rot90(m,2)), D:-1:1);
end