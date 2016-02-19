classdef Curvilinear < grid.Structured
    % General grid mapping
    methods (Abstract)
        % baseGrid returns the domain grid of the mapping.
        g = baseGrid(obj);

        % mapping returns the mapped coordinates as a grid.Function
        m = mapping(obj);
    end
end