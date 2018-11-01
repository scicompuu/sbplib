classdef Nodes < grid.Grid
    properties
        coords
    end

    methods
        % Creates a grid with one point for each row in coords.
        % The dimension equals the number of columns in coords.
        function obj = Nodes(coords)
            obj.coords = coords;
        end

        function o = N(obj)
            o = size(obj.coords, 1);
        end

        % d returns the spatial dimension of the grid
        function o = D(obj)
            o = size(obj.coords, 2);
        end

        % points returns a n x d matrix containing the coordinates for all points.
        function X = points(obj)
            X = obj.coords;
        end

        % Restricts the grid function gf on obj to the subgrid g.
        function gf = restrictFunc(obj, gf, g)
            error('Not implemented');
        end

        % Projects the grid function gf on obj to the grid g.
        function gf = projectFunc(obj, gf, g)
            error('Not implemented');
        end

        % Return the grid.boundaryIdentifiers of all boundaries in a cell array.
        function bs = getBoundaryNames(obj)
            error('Not implemented');
        end

        % Return coordinates for the given boundary
        function b = getBoundary(obj, name)
            error('Not implemented');
        end
    end
end
