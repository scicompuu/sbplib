classdef Grid < grid.Grid
    properties
        grids
        connections
        boundaryGroups

        nPoints
    end

    % General multiblock grid
    methods

        % grids -- cell array of N grids
        % connections -- NxN upper triangular cell matrix. connections{i,j}
        %                specifies the connection between block i and j. If
        %                it's empty there is no connection otherwise it's a 2
        %                -cell-vector with strings naming the boundaries to be
        %                connected. (inverted coupling?)
        %% Should we have boundary groups at all? maybe it can be handled in a
        %% cleaner way outside of the class. Maybe the grid class doesn't need to care at all. All boundaryGroup interaction is in DiffOp?
        function obj = Grid(grids, connections, boundaryGroups)
            obj.grids = grids;
            obj.connections = connections;

            obj.nPoints = 0;
            for i = 1:length(grids)
                obj.nPoints = obj.nPoints + grids{i}.N();
            end

            % if iscell(boundaryGroups)
        end

        function n = size(obj)
            n = length(obj.grids);
        end

        % N returns the number of points in the grid
        function o = N(obj)
            o = obj.nPoints;
        end

        % Ns returns the number of points in each sub grid as a vector
        function o = Ns(obj)
            ns = zeros(1,obj.nBlocks);
            for i = 1:obj.nBlocks;
                ns(i) = obj.grids{i}.N();
            end
            o = ns;
        end

        function n = nBlocks(obj)
            n = length(obj.grids);
        end

        % d returns the spatial dimension of the grid
        function o = D(obj)
            o = obj.grids{1}.D();
        end

        % points returns a n x d matrix containing the coordinates for all points.
        function X = points(obj)
            X = [];
            for i = 1:length(obj.grids)
                X = [X; obj.grids{i}.points];
            end
        end

        % Split a grid function on obj to a cell array of grid function on each block
        function gfs = splitFunc(obj, gf)
            nComponents = length(gf)/obj.nPoints;
            nBlocks = length(obj.grids);

            % Collect number of points in each block
            N = zeros(1,nBlocks);
            for i = 1:nBlocks
                N(i) = obj.grids{i}.N();
            end

            gfs = mat2cell(gf, N, 1);
        end

        % Converts a gridfunction to a set of plot matrices
        % Takes a grid function and and a structured grid.
        function F = funcToPlotMatrices(obj, gf)
            gfs = obj.splitFunc(gf);

            F = cell(1, obj.nBlocks());

            for i = 1:obj.nBlocks()
                F{i} = grid.funcToPlotMatrix(obj.grids{i}, gfs{i});
            end
        end


        % Restricts the grid function gf on obj to the subgrid g.
        function gf = restrictFunc(obj, gf, g)
            gfs = obj.splitFunc(gf);

            for i = 1:length(obj.grids)
                gfs{i} = obj.grids{i}.restrictFunc(gfs{i}, g.grids{i});
            end

            gf = cell2mat(gfs);
        end

        % Projects the grid function gf on obj to the grid g.
        function o = projectFunc(obj, gf, g)
            error('not implemented')

            p = g.points();
            o = zeros(length(p),1);
            for i = 1:length(p)
                I = whatGrid(p(i));
                o(i) = obj.grids{I}.projectFunc(gf, p(i));
            end


            function I = whatGrid(p)
                % Find what grid a point lies on
            end

        end

        function bs = getBoundaryNames(obj)
            bs = [];
        end

        % Return coordinates for the given boundary
        function b = getBoundary(obj, name)
            b = [];
        end
    end
end
