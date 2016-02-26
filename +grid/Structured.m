classdef Structured < grid.Grid
    methods (Abstract)
        % Returns the size of the grid in each dimension m = [mx my mz ...]
        m = size(obj);
        h = scaling(obj);
    end
end
