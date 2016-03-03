% Converts a gridfunction to a plot matrix
% Takes a grid function and and a structured grid.
function F = funcToPlotMatrix(g, gf)
    reshapeToPlotMatrix(gf, g.size());
end