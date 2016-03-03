classdef Cartesian < grid.Structured
    properties
        n % Number of points in the grid
        d % Number of dimensions
        m % Number of points in each direction
        x % Cell array of vectors with node placement for each dimension.
        h % Spacing/Scaling
    end

    % General d dimensional grid with n points
    methods
        % Creates a cartesian grid given vectors conatining the coordinates
        % in each direction
        function obj = Cartesian(varargin)
            obj.d = length(varargin);

            for i = 1:obj.d
                obj.x{i} = varargin{i};
                obj.m(i) = length(varargin{i});
            end

            obj.n = prod(obj.m);
            if obj.n == 0
                error('grid:Cartesian:EmptyGrid','Input parameter gives an empty grid.')
            end

            obj.h = [];
        end
        % n returns the number of points in the grid
        function o = N(obj)
            o = obj.n;
        end

        % d returns the spatial dimension of the grid
        function o = D(obj)
            o = obj.d;
        end

        function m = size(obj)
            m = obj.m;
        end

        % points returns a n x d matrix containing the coordianets for all points.
        % points are ordered according to the kronecker product with X*Y*Z
        function X = points(obj)
            X = zeros(obj.n, obj.d);

            for i = 1:obj.d
                if iscolumn(obj.x{i})
                    c = obj.x{i};
                else
                    c = obj.x{i}';
                end

                m_before = prod(obj.m(1:i-1));
                m_after = prod(obj.m(i+1:end));

                X(:,i) = kr(ones(m_before,1),c,ones(m_after,1));
            end
        end

        % matrices returns a cell array with coordinates in matrix form.
        % For 2d case these will have to be transposed to work with plotting routines.
        function X = matrices(obj)

            if obj.d == 1 % There is no 1d matrix data type in matlab, handle special case
                X{1} = reshape(obj.x{1}, [obj.m 1]);
                return
            end

            X = cell(1,obj.d);
            for i = 1:obj.d
                s = ones(1,obj.d);
                s(i) = obj.m(i);

                t = reshape(obj.x{i},s);

                s = obj.m;
                s(i) = 1;
                X{i} = repmat(t,s);
            end
        end

        function h = scaling(obj)
            if isempty(obj.h)
                error('grid:Cartesian:NoScalingSet', 'No scaling set')
            end

            h = obj.h;
        end

        % Restricts the grid function gf on obj to the subgrid g.
        % Only works for even multiples
        function gf = restrictFunc(obj, gf, g)
            m1 = obj.m;
            m2 = g.m;

            % Check the input
            if prod(m1) ~= numel(gf)
                error('grid:Cartesian:restrictFunc:NonMatchingFunctionSize', 'The grid function has to few or too many points.');
            end

            if ~all(mod(m1-1,m2-1) == 0)
                error('grid:Cartesian:restrictFunc:NonMatchingGrids', 'Only integer downsamplings are allowed');
            end

            % Calculate stride for each dimension
            stride = (m1-1)./(m2-1);

            % Create downsampling indecies
            I = {};
            for i = 1:length(m1)
                I{i} = 1:stride(i):m1(i);
            end

            gf = reshapeRowMaj(gf, m1);
            gf = gf(I{:});
            gf = reshapeRowMaj(gf, prod(m2));
        end

        % Projects the grid function gf on obj to the grid g.
        function gf = projectFunc(obj, gf, g)
            error('grid:Cartesian:NotImplemented')
        end
    end
end