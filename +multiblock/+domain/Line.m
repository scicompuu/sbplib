classdef Line < multiblock.Definition
    properties

    xlims
    blockNames % Cell array of block labels
    nBlocks
    connections % Cell array specifying connections between blocks
    boundaryGroups % Structure of boundaryGroups

    end


    methods
        % Creates a divided line
        % x is a vector of boundary and interface positions.
        % blockNames: cell array of labels. The id is default.
        function obj = Line(x,blockNames)
            default_arg('blockNames',[]);

            N = length(x)-1; % number of blocks in the x direction.

            if ~issorted(x)
                error('The elements of x seem to be in the wrong order');
            end

            % Dimensions of blocks and number of points
            blockTi = cell(N,1);
            xlims = cell(N,1);
            for i = 1:N
                xlims{i} = {x(i), x(i+1)};
            end

            % Interface couplings
            conn = cell(N,N);
            for i = 1:N
                conn{i,i+1} = {'r','l'};
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
            L = { {1, 'l'} };
            R = { {N, 'r'} };
            boundaryGroups.L = multiblock.BoundaryGroup(L);
            boundaryGroups.R = multiblock.BoundaryGroup(R);
            boundaryGroups.all = multiblock.BoundaryGroup([L,R]);

            obj.connections = conn;
            obj.nBlocks = nBlocks;
            obj.boundaryGroups = boundaryGroups;
            obj.xlims = xlims;

        end


        % Returns a multiblock.Grid given some parameters
        % ms: cell array of m values 
        % For same m in every block, just input one scalar.
        function g = getGrid(obj, ms, varargin)

            default_arg('ms',21)

            % Extend ms if input is a single scalar
            if (numel(ms) == 1) && ~iscell(ms)
                m = ms;
                ms = cell(1,obj.nBlocks);
                for i = 1:obj.nBlocks
                    ms{i} = m;
                end
            end

            grids = cell(1, obj.nBlocks);
            for i = 1:obj.nBlocks
                grids{i} = grid.equidistant(ms{i}, obj.xlims{i}, obj.ylims{i});
            end

            g = multiblock.Grid(grids, obj.connections, obj.boundaryGroups);
        end

        % Returns a multiblock.Grid given some parameters
        % ms: cell array of m values 
        % For same m in every block, just input one scalar.
        function g = getStaggeredGrid(obj, ms, varargin)

            default_arg('ms',21)

            % Extend ms if input is a single scalar
            if (numel(ms) == 1) && ~iscell(ms)
                m = ms;
                ms = cell(1,obj.nBlocks);
                for i = 1:obj.nBlocks
                    ms{i} = m;
                end
            end

            grids = cell(1, obj.nBlocks);
            for i = 1:obj.nBlocks
                [g_primal, g_dual] = grid.primalDual1D(ms{i}, obj.xlims{i});
                grids{i} = grid.Staggered1d(g_primal, g_dual);
            end

            g = multiblock.Grid(grids, obj.connections, obj.boundaryGroups);
        end

        % label is the type of label used for plotting,
        % default is block name, 'id' show the index for each block.
        function show(obj, label)
            default_arg('label', 'name')

            m = 10;
            figure
            for i = 1:obj.nBlocks
               x = linspace(obj.xlims{i}{1}, obj.xlims{i}{2}, m);
               y = 0*x + 0.05* ( (-1)^i + 1 ) ;
               plot(x,y,'+');
               hold on 
            end
            hold off

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

            legend(labels)
            axis equal
        end

        % Returns the grid size of each block in a cell array
        % The input parameters are determined by the subclass
        function ms = getGridSizes(obj, varargin)
        end
    end
end
