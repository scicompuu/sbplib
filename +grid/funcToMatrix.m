% Converts a gridfunction to a matrix
% Takes a grid function and and a structured grid.
function F = funcToMatrix(g, gf)
    reshapeRowMaj(gf, g.size());
end