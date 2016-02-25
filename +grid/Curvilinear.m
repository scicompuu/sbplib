classdef Curvilinear < grid.Structured & grid.Mapped
    properties
        logic % Grid of Logical domain
        coords % N x D matrix with coordinates of each point in the physical domain
    end

    methods
        % Creates a curvilinear grid.
        % Ex: grid.Curvilinear(mapping, xi, eta, ...)
        %    mapping     -- either a matrix or a cell array with physical coordinates.
        %                   A matrix should be a grid function (N*D x 1 vector) or a N x D
        %                   A cell array should be a 1 x D cell array with either N x 1 vectors
        %                   or matrices of the same dimesions as the logical grid.
        %   xi, eta, ... -- are the coordinate positions of the cartesian logical grid.
        function obj = Curvilinear(mapping, varargin)
            xi = varargin;
            obj.logic = grid.Cartesian(xi{:});

            % If mapping is a function evaluate it
            if isa(mapping, 'function_handle')
                mapping = grid.evalOn(obj.logic, mapping);
            end

            D = obj.logic.D();
            N = obj.logic.N();

            obj.coords = zeros(N,D);

            if iscell(mapping)
                if ~isequal(size(mapping),[1 D])
                    error('grid:Curvilinear:Curvilinear','The cell array must be a row array.');
                end

                if isequal(size(mapping{1}),[N 1])
                    obj.coords = cell2mat(mapping);
                elseif isequal(size(mapping{1}), obj.logic.m)
                    for i = 1:length(mapping)
                        obj.coords(:,i) = reshapeRowMaj(mapping{i}, [N 1]);
                    end
                else
                    error('grid:Curvilinear:Curvilinear','The matrix must have size [N 1] or the same dimension as the grid. Actual: %s', toString(obj.logic.m));
                end

            elseif isnumeric(mapping)
                if isequal(size(mapping), [N, D])
                    obj.coords = mapping;
                elseif isequal(size(mapping), [N*D, 1])
                    obj.coords = reshapeRowMaj(mapping,[N D]);
                else
                    error('grid:Curvilinear:Curvilinear','A matrix mapping must be of size [N D] or [N*D 1].');
                end
            else
                error('grid:Curvilinear:Curvilinear','mapping must be a matrix or a cell array.');
            end
        end

        function m = size(obj)
            m = obj.logic.size();
        end

        % logicalGrid returns the domain grid of the mapping.
        function g = logicalGrid(obj)
            g = obj.logic;
        end

        % mapping returns the mapped coordinates as a grid.Function
        function m = mapping(obj);
            m = obj.coords;
        end

        % n returns the number of points in the grid
        function o = N(obj)
            o = obj.logic.N();
        end

        % d returns the spatial dimension of the grid
        function o = D(obj)
            o = obj.logic.D();
        end

        % points returns a n x d matrix containing the coordinates for all points.
        function X = points(obj)
            X = obj.coords;
        end

        % Restricts the grid function gf on obj to the subgrid g.
        function gf = restrictFunc(obj, gf, g)
            gf = obj.logic.restrictFunc(gf, g.baseGrid());
        end

        % Projects the grid function gf on obj to the grid g.
        function gf = projectFunc(obj, gf, g)
            gf = obj.logic.projectFunc(gf,g.baseGrid());
        end
    end
end