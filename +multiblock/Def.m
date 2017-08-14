classdef Def
    properties
        nBlocks
        blockMaps % Maps from logical blocks to physical blocks build from transfinite interpolation
        blockNames
        connections % Cell array specifying connections between blocks
        boundaryGroups % Structure of boundaryGroups
    end

    methods
        % Defines a multiblock setup for transfinite interpolation blocks
        % TODO: How to bring in plotting of points?
        function obj = Def(blockMaps, connections, boundaryGroups, blockNames)
            default_arg('boundaryGroups', struct());
            default_arg('blockNames',{});

            nBlocks = length(blockMaps);

            obj.nBlocks = nBlocks;

            obj.blockMaps = blockMaps;

            assert(all(size(connections) == [nBlocks, nBlocks]));
            obj.connections = connections;


            if isempty(blockNames)
                obj.blockNames = cell(1, nBlocks);
                for i = 1:length(blockMaps)
                    obj.blockNames{i} = sprintf('%d', i);
                end
            else
                assert(length(blockNames) == nBlocks);
                obj.blockNames = blockNames;
            end

            obj.boundaryGroups = boundaryGroups;
        end

        function g = getGrid(obj, varargin)
            ms = obj.getGridSizes(varargin{:});

            grids = cell(1, obj.nBlocks);
            for i = 1:obj.nBlocks
                grids{i} = grid.equidistantCurvilinear(obj.blockMaps{i}.S, ms{i});
            end

            g = multiblock.Grid(grids, obj.connections, obj.boundaryGroups);
        end

        function show(obj, label, gridLines, varargin)
            default_arg('label', true)
            default_arg('gridLines', false);

            if ~label && ~gridLines
                for i = 1:obj.nBlocks
                    obj.blockMaps{i}.show(2,2);
                end
                axis equal
                return
            end

            if gridLines
                ms = obj.getGridSizes(varargin{:});
                for i = 1:obj.nBlocks
                    obj.blockMaps{i}.show(ms{i}(1),ms{i}(2));
                end
            end

            if label
                for i = 1:obj.nBlocks
                    parametrization.Ti.label(obj.blockMaps{i}, obj.blockNames{i});
                end
            end

            axis equal
        end
    end

    methods (Abstract)
        % Returns the grid size of each block in a cell array
        % The input parameters are determined by the subclass
        ms = getGridSizes(obj, varargin)
        % end
    end

end


