classdef Definition
    methods (Abstract)

        % Returns a multiblock.Grid given some parameters
        g = getGrid(obj, varargin)

        % label is the type of label used for plotting,
        % default is block name, 'id' show the index for each block.
        show(obj, label, gridLines, varargin)

        % Returns the grid size of each block in a cell array
        % The input parameters are determined by the subclass
        ms = getGridSizes(obj, varargin)
    end
end
