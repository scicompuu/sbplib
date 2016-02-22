% Takes a grid function and reshapes it into a matrix of shape m for plotting.
function F = reshapeToPlotMatrix(gf, m)
    D = length(m);

    switch D
        case 1
            F = gf;
        case 2
            F = reshape(gf, rot90(m,2));
        case 3
            % After the reshape the indecies will be M(z,y,x). Plot need them to be M(y,x,z)
            p = [2 3 1]; % Permuation
            F = permute(reshape(gf,rot90(m,2)), p);
        otherwise
            error('reshapeToPlotMatrix:NotImplemented','Grid function to matrix is not implemented for dimension = %d', length(m));
    end
end