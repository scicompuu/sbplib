classdef Grid < grid.Grid
    % General multiblock grid
    methods (Abstract)
        % NBlocks returns the number of blocks in the grid.
        o = NBlocks(obj);

        % Grid returns the ith grid in the multiblockgrid
        gs = Grid(obj,i);

        % Grids returns a cell array of all the grids in the multiblock grid.
        gs = Grids(obj);
    end
end