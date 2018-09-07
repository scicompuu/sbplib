classdef Mapped < grid.Grid
    % General grid mapping
    methods (Abstract)
        % logicalGrid returns the domain grid of the mapping.
        g = logicalGrid(obj);

        % mapping returns the mapped coordinates as a N x D component matrix
        m = mapping(obj);
    end
end