classdef Grid < handle
    % General d dimensional grid with n points
    methods (Abstract)
        % n returns the number of points in the grid
        o = N(obj)

        % d returns the spatial dimension of the grid
        o = D(obj)

        % points returns a n x d matrix containing the coordinates for all points.
        X = points(obj)

        % Restricts the grid function gf on obj to the subgrid g.
        gf = restrictFunc(obj, gf, g)

        % Projects the grid function gf on obj to the grid g.
        gf = projectFunc(obj, gf, g)

        % Return the names of all boundaries in this grid.
        bs = getBoundaryNames(obj)

        % Return coordinates for the given boundary
        b = getBoundary(obj, name)
    end
end
