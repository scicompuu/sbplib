classdef Rectangle < multiblock.Definition
    properties

    end


    methods
        % Creates a ...
        % x and y are vectors of boundary and interface positions.
        function obj = Rectangle(x,y)

        end


        % Returns a multiblock.Grid given some parameters
        function g = getGrid(obj, varargin)
        end

        % label is the type of label used for plotting,
        % default is block name, 'id' show the index for each block.
        function show(obj, label, gridLines, varargin)

        end

        % Returns the grid size of each block in a cell array
        % The input parameters are determined by the subclass
        function ms = getGridSizes(obj, varargin)
        end
    end
end
