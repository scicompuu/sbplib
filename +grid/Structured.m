classdef Structured < grid.Grid
    methods (Abstract)
        % Returns the size of the grid in each dimension m = [mx my mz ...]
        m = size(obj); % Is this a good idea? Isn't immersed a structured grid?
        h = scaling(obj);
    end
end
