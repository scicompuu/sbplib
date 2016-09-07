classdef Curvilinear < grid.Structured & grid.Mapped
    properties
        logic % Grid of Logical domain
        coords % N x D matrix with coordinates of each point in the physical domain
    end

    methods
        % Creates a curvilinear grid.
        % Ex: grid.Curvilinear(mapping, xi, eta, ...)
        %    mapping     -- either a function handle, a matrix or a cell array with physical coordinates.
        %                   A function handle should be a vector valued function of the coordinate mapping.
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
                obj.coords = cellMappingToCoords(mapping, N, D, obj.logic.m);
            elseif isnumeric(mapping)
                obj.coords = matrixMappingToCoords(mapping, N, D);
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

        function h = scaling(obj)
            if isempty(obj.logic.h)
                error('grid:Curvilinear:NoScalingSet','No scaling set');
            end
            h = obj.logic.h;
        end

        % Return the names of all boundaries in this grid.
        function bs = getBoundaryNames(obj)
            bs = obj.logic.getBoundaryNames();
        end

        % Return coordinates for the given boundary
        function X = getBoundary(obj, name)
              % In what dimension is the boundary?
            switch name
                case {'l', 'r', 'w', 'e'}
                    D = 1;
                case {'s', 'n'}
                    D = 2;
                case {'d', 'u'}
                    D = 3;
                otherwise
                    error('not implemented');
            end

            % At what index is the boundary?
            switch name
                case {'l', 'w', 's', 'd'}
                    index = 1;
                case {'r', 'e', 'n', 'u'}
                    index = obj.logic.m(D);
                otherwise
                    error('not implemented');
            end



            I = cell(1, obj.D);
            for i = 1:obj.D
                if i == D
                    I{i} = index;
                else
                    I{i} = ':';
                end
            end

            % Calculate size of result:
            m = obj.logic.m;
            m(D) = [];
            N = prod(m);

            X = zeros(N, obj.D);

            p = obj.points;
            for i = 1:obj.D()
                coordMat{i} = reshapeRowMaj(p(:,i), obj.logic.m);
            end

            for i = 1:length(coordMat)
                Xtemp = coordMat{i}(I{:});
                X(:,i) = reshapeRowMaj(Xtemp, [N,1]);
            end
        end
    end
end


function coords = cellMappingToCoords(mapping, N, D, m)
    if ~isequal(size(mapping),[1 D])
        error('grid:Curvilinear:Curvilinear','The cell array must be a 1xD array.');
    end

    if isequal(size(mapping{1}),[N 1])
        coords = cell2mat(mapping);
    elseif isequal(size(mapping{1}), m)
        for i = 1:length(mapping)
            coords(:,i) = reshapeRowMaj(mapping{i}, [N 1]);
        end
    else
        error('grid:Curvilinear:Curvilinear','The matrix must have size [N 1] or the same dimension as the grid. Actual: %s', toString(m));
    end
end

function coords = matrixMappingToCoords(mapping, N, D)
    if isequal(size(mapping), [N, D])
        coords = mapping;
    elseif isequal(size(mapping), [N*D, 1])
        coords = reshapeRowMaj(mapping,[N D]);
    else
        error('grid:Curvilinear:Curvilinear','A matrix mapping must be of size [N D] or [N*D 1].');
    end
end
