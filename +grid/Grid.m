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
        gf = restrictFunc(gf, g)

        % Projects the grid function gf on obj to the grid g.
        gf = projectFunc(gf, g)
    end
end


%% Should it be able to return a cell size aswell? For an equidistant grid this would be know
%% for other grids the constructor would have to make something up.
%% For example the grid.Cartesian constructor would take a h (1 x d) vector as an in parameter.


%Should define boundaries somehow, look in stitchSchemes.