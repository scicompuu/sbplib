classdef Rectangle < multiblock.Definition
    properties

    blockTi % Transfinite interpolation objects used for plotting
    xlims
    ylims
    blockNames % Cell array of block labels
    nBlocks
    connections % Cell array specifying connections between blocks
    boundaryGroups % Structure of boundaryGroups

    end


    methods
        % Creates a divided rectangle
        % x and y are vectors of boundary and interface positions.
        % blockNames: cell array of labels. The id is default.
        function obj = Rectangle(x,y,blockNames)
            default_arg('blockNames',[]);

            n = length(y)-1; % number of blocks in the y direction.
            m = length(x)-1; % number of blocks in the x direction.
            N = n*m; % number of blocks

            if ~issorted(x)
                error('The elements of x seem to be in the wrong order');
            end
            if ~issorted(flip(y))
                error('The elements of y seem to be in the wrong order');
            end

            % Dimensions of blocks and number of points
            blockTi = cell(N,1);
            xlims = cell(N,1);
            ylims = cell(N,1);
            for i = 1:n
                for j = 1:m
                    p1 = [x(j), y(i+1)];
                    p2 = [x(j+1), y(i)];
                    I = flat_index(m,j,i);
                    blockTi{I} = parametrization.Ti.rectangle(p1,p2);
                    xlims{I} = {x(j), x(j+1)};
                    ylims{I} = {y(i+1), y(i)};
                end
            end

            % Interface couplings
            conn = cell(N,N);
            for i = 1:n
                for j = 1:m
                    I = flat_index(m,j,i);
                    if i < n
                        J = flat_index(m,j,i+1);
                        conn{I,J} = {'s','n'};
                    end

                    if j < m
                        J = flat_index(m,j+1,i);
                        conn{I,J} = {'e','w'};
                    end
                end
            end

            % Block names (id number as default)
            if isempty(blockNames)
                obj.blockNames = cell(1, N);
                for i = 1:N
                    obj.blockNames{i} = sprintf('%d', i);
                end
            else
                assert(length(blockNames) == N);
                obj.blockNames = blockNames;
            end
            nBlocks = N;

            % Boundary groups
            boundaryGroups = struct();
            nx = m;
            ny = n;
            E = cell(1,ny);
            W = cell(1,ny);
            S = cell(1,nx);
            N = cell(1,nx);
            for i = 1:ny
                E_id = flat_index(m,nx,i);
                W_id = flat_index(m,1,i);
                E{i} = {E_id,'e'};
                W{i} = {W_id,'w'};
            end
            for j = 1:nx
                S_id = flat_index(m,j,ny);
                N_id = flat_index(m,j,1);
                S{j} = {S_id,'s'};
                N{j} = {N_id,'n'};
            end
            boundaryGroups.E = multiblock.BoundaryGroup(E);
            boundaryGroups.W = multiblock.BoundaryGroup(W);
            boundaryGroups.S = multiblock.BoundaryGroup(S);
            boundaryGroups.N = multiblock.BoundaryGroup(N);
            boundaryGroups.all = multiblock.BoundaryGroup([E,W,S,N]);
            boundaryGroups.WS = multiblock.BoundaryGroup([W,S]);
            boundaryGroups.WN = multiblock.BoundaryGroup([W,N]);
            boundaryGroups.ES = multiblock.BoundaryGroup([E,S]);
            boundaryGroups.EN = multiblock.BoundaryGroup([E,N]);

            obj.connections = conn;
            obj.nBlocks = nBlocks;
            obj.boundaryGroups = boundaryGroups;
            obj.blockTi = blockTi;
            obj.xlims = xlims;
            obj.ylims = ylims;

        end


        % Returns a multiblock.Grid given some parameters
        % ms: cell array of [mx, my] vectors
        % For same [mx, my] in every block, just input one vector.
        % Currently defaults to equidistant grid if varargin is empty.
        % If varargin is non-empty, the first argument should supply the grid type, followed by
        % additional arguments required to construct the grid.
        % Grid types:
        %          'equidist' - equidistant grid
        %                       Additional argumets: none
        %          'boundaryopt' - boundary optimized grid based on boundary
        %                          optimized SBP operators
        %                          Additional arguments: order, stencil option
        % Example: g = getGrid() - the local blocks are 21x21 equidistant grids.
        %          g = getGrid(ms,) - block i is an equidistant grid with size given by ms{i}.
        %          g = getGrid(ms,'equidist') - block i is an equidistant grid with size given by ms{i}.
        %          g = getGrid(ms,'boundaryopt',4,'minimal') - block i is a Cartesian grid with size given by ms{i}
        %              and nodes placed according to the boundary optimized minimal 4th order SBP operator.
        function g = getGrid(obj, ms, varargin)

            default_arg('ms',[21,21])

            % Extend ms if input is a single vector
            if (numel(ms) == 2) && ~iscell(ms)
                m = ms;
                ms = cell(1,obj.nBlocks);
                for i = 1:obj.nBlocks
                    ms{i} = m;
                end
            end
            if isempty(varargin) || strcmp(varargin{1},'equidist')
               gridgenerator = @(m,xlim,ylim)grid.equidistant(m,xlim,ylim);
            elseif strcmp(varargin{1},'boundaryopt')
                order = varargin{2};
                stenciloption = varargin{3};
                gridgenerator = @(m,xlim,ylim)grid.boundaryOptimized(m,xlim,ylim,...
                    order,stenciloption);
            else
                error('No grid type supplied!');
            end
            grids = cell(1, obj.nBlocks);
            for i = 1:obj.nBlocks
                grids{i} = gridgenerator(ms{i}, obj.xlims{i}, obj.ylims{i});
            end

            g = multiblock.Grid(grids, obj.connections, obj.boundaryGroups);
        end

        % label is the type of label used for plotting,
        % default is block name, 'id' show the index for each block.
        function show(obj, label, gridLines, varargin)
            default_arg('label', 'name')
            default_arg('gridLines', false);

            if isempty('label') && ~gridLines
                for i = 1:obj.nBlocks
                    obj.blockTi{i}.show(2,2);
                end
                axis equal
                return
            end

            if gridLines
                m = 10;
                for i = 1:obj.nBlocks
                    obj.blockTi{i}.show(m,m);
                end
            end


            switch label
                case 'name'
                    labels = obj.blockNames;
                case 'id'
                    labels = {};
                    for i = 1:obj.nBlocks
                        labels{i} = num2str(i);
                    end
                otherwise
                    axis equal
                    return
            end

            for i = 1:obj.nBlocks
                parametrization.Ti.label(obj.blockTi{i}, labels{i});
            end

            axis equal
        end

        % Returns the grid size of each block in a cell array
        % The input parameters are determined by the subclass
        function ms = getGridSizes(obj, varargin)
        end
    end
end
