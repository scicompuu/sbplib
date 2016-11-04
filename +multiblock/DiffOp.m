classdef DiffOp < scheme.Scheme
    properties
        grid
        order
        diffOps
        D
        H

        blockmatrixDiv
    end

    methods
        function obj = DiffOp(doHand, grid, order, doParam)
            %  doHand -- may either be a function handle or a cell array of
            %            function handles for each grid. The function handle(s)
            %            should be on the form do = doHand(grid, order, ...)
            %            Additional parameters for each doHand may be provided in
            %            the doParam input.
            %    grid -- a multiblock grid
            %   order -- integer specifying the order of accuracy
            % doParam -- may either be a cell array or a cell array of cell arrays
            %            for each block. If it is a cell array with length equal
            %            to the number of blocks then each element is sent to the
            %            corresponding function handle as extra parameters:
            %            doHand(..., doParam{i}{:}) Otherwise doParam is sent as
            %            extra parameters to all doHand: doHand(..., doParam{:})
            default_arg('doParam', [])

            [getHand, getParam] = parseInput(doHand, grid, doParam);

            nBlocks = grid.nBlocks();

            obj.order = order;

            % Create the diffOps for each block
            obj.diffOps = cell(1, nBlocks);
            for i = 1:nBlocks
                h = getHand(i);
                p = getParam(i);
                if ~iscell(p)
                    p = {p};
                end
                obj.diffOps{i} = h(grid.grids{i}, order, p{:});
            end


            % Build the norm matrix
            H = cell(nBlocks, nBlocks);
            for i = 1:nBlocks
                H{i,i} = obj.diffOps{i}.H;
            end
            obj.H = blockmatrix.toMatrix(H);


            % Build the differentiation matrix
            obj.blockmatrixDiv = {grid.Ns, grid.Ns};
            D = blockmatrix.zero(obj.blockmatrixDiv);
            for i = 1:nBlocks
                D{i,i} = obj.diffOps{i}.D;
            end

            for i = 1:nBlocks
                for j = 1:nBlocks
                    intf = grid.connections{i,j};
                    if isempty(intf)
                        continue
                    end


                    [ii, ij] = obj.diffOps{i}.interface(intf{1}, obj.diffOps{j}, intf{2});
                    D{i,i} = D{i,i} + ii;
                    D{i,j} = D{i,j} + ij;

                    [jj, ji] = obj.diffOps{j}.interface(intf{2}, obj.diffOps{i}, intf{1});
                    D{j,j} = D{j,j} + jj;
                    D{j,i} = D{j,i} + ji;
                end
            end
            obj.D = blockmatrix.toMatrix(D);


            function [getHand, getParam] = parseInput(doHand, grid, doParam)
                if ~isa(grid, 'multiblock.Grid')
                    error('multiblock:DiffOp:DiffOp:InvalidGrid', 'Requires a multiblock grid.');
                end

                if iscell(doHand) && length(doHand) == grid.nBlocks()
                    getHand = @(i)doHand{i};
                elseif isa(doHand, 'function_handle')
                    getHand = @(i)doHand;
                else
                    error('multiblock:DiffOp:DiffOp:InvalidGridDoHand', 'doHand must be a function handle or a cell array of length grid.nBlocks');
                end

                if isempty(doParam)
                    getParam = @(i){};
                elseif iscell(doParam) && length(doParam) == grid.nBlocks()
                    getParam = @(i)doParam{i};
                else
                    getParam = @(i)doParam;
                end
            end
        end

        function ops = splitOp(obj, op)
            % Splits a matrix operator into a cell-matrix of matrix operators for
            % each grid.
            ops = sparse2cell(op, obj.NNN);
        end

        % Creates the closure and penalty matrix for a given boundary condition,
        %    boundary -- the name of the boundary on the form [id][name] where
        %                id is the number of a block and name is the name of a
        %                boundary of that block example: 1s or 3w
        function [closure, penalty] = boundary_condition(obj, boundary, type)
            I = boundary{1};
            name = boundary{2};

            % Get the closure and penaly matrices
            [blockClosure, blockPenalty] = obj.diffOps{I}.boundary_condition(name, type);

            % Expand to matrix for full domain.
            div = obj.blockmatrixDiv;
            closure = blockmatrix.zero(div);
            closure{I,I} = blockClosure;
            closure = blockmatrix.toMatrix(closure);

            div{2} = size(blockPenalty, 2); % Penalty is a column vector
            if ~iscell(blockPenalty)
                p = blockmatrix.zero(div);
                p{I} = blockPenalty;
                penalty = blockmatrix.toMatrix(p);
            else
                for i = 1:length(blockPenalty)
                    p = blockmatrix.zero(div);
                    p{I} = blockPenalty{i};
                    penalty{i} = blockmatrix.toMatrix(p);
                end
            end
        end

        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary)

        end

        % Size returns the number of degrees of freedom
        function N = size(obj)
            N = 0;
            for i = 1:length(obj.diffOps)
                N = N + obj.diffOps{i}.size();
            end
        end
    end
end
