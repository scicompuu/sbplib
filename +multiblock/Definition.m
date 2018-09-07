classdef Definition
    methods (Abstract)

        % Returns a multiblock.Grid given some parameters
        g = getGrid(obj, varargin)

        % label is the type of label used for plotting,
        % default is block name, 'id' show the index for each block.
        show(obj, label, gridLines, varargin)
    end
end
