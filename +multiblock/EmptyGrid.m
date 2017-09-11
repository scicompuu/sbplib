classdef EmptyGrid < grid.Empty & multiblock.Grid
    methods
        function obj = EmptyGrid(D)
            obj@multiblock.Grid({},cell(0,0));
            obj@grid.Empty(D);
        end

        % n returns the number of points in the grid
        function o = N(obj)
            o = N@grid.Empty(obj);
        end

        % d returns the spatial dimension of the grid
        function o = D(obj)
            o = D@grid.Empty(obj);
        end

        % points returns a n x d matrix containing the coordinates for all points.
        function X = points(obj)
            X = points@grid.Empty(obj);
        end

        % Restricts the grid function gf on obj to the subgrid g.
        function gf = restrictFunc(obj, gf, g)
            gf = restrictFunc@grid.Empty(obj);
        end

        % Projects the grid function gf on obj to the grid g.
        function gf = projectFunc(obj, gf, g)
            gf = projectFunc@grid.Empty(obj);
        end

        % Return the grid.boundaryIdentifiers of all boundaries in a cell array.
        function bs = getBoundaryNames(obj)
            bs = getBoundaryNames@grid.Empty(obj);
        end

        % Return coordinates for the given boundary
        function b = getBoundary(obj, name)
            b = getBoundary@grid.Empty(name);
        end

        function s = size(obj)
            s = size@multiblock.Grid(obj);
        end
    end
end
