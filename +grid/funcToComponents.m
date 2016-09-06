% funcToComponents converts a grid function to a N x k matrix, where
% k is the number of vector components of the gridfunction and N is the
% number of points in the grid.
%
% Takes a grid function and and a grid.
function F = funcToComponents(g, gf);
    F = reshapeRowMaj(gf, [g.N, length(gf)/g.N]);
end