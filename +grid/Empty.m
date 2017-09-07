classdef Empty < grid.Grid & grid.Structured
    properties
        dim
    end

    methods
        function obj = Empty(D)
            obj.dim = D;
        end
        % n returns the number of points in the grid
        function o = N(obj)
            o = 0;
        end

        % d returns the spatial dimension of the grid
        function o = D(obj)
            o = obj.dim;
        end

        % points returns a n x d matrix containing the coordinates for all points.
        function X = points(obj)
            X = sparse(0,obj.dim);
        end

        % Restricts the grid function gf on obj to the subgrid g.
        function gf = restrictFunc(obj, gf, g)
            error('Restrict does not make sense for an empty grid')
        end

        % Projects the grid function gf on obj to the grid g.
        function gf = projectFunc(obj, gf, g)
            error('Project does not make sense for an empty grid')
        end

        % Return the grid.boundaryIdentifiers of all boundaries in a cell array.
        function bs = getBoundaryNames(obj)
            bs = {};
        end

        % Return coordinates for the given boundary
        function b = getBoundary(obj, name)
            b = sparse(0,obj.dim-1);
        end

        function h = scaling(obj)
            h = 1;
        end

        function s = size(obj)
            s = zeros(1, obj.dim);
        end
    end
end