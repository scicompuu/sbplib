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
        function obj = DiffOp(doHand, g, order, doParam)
            %  doHand -- may either be a function handle or a cell array of
            %            function handles for each grid. The function handle(s)
            %            should be on the form do = doHand(grid, order, ...)
            %            Additional parameters for each doHand may be provided in
            %            the doParam input.
            %       g -- a multiblock grid
            %   order -- integer specifying the order of accuracy
            % doParam -- may either be a cell array or a cell array of cell arrays
            %            for each block. If it is a cell array with length equal
            %            to the number of blocks then each element is sent to the
            %            corresponding function handle as extra parameters:
            %            doHand(..., doParam{i}{:}) Otherwise doParam is sent as
            %            extra parameters to all doHand: doHand(..., doParam{:})
            default_arg('doParam', [])

            [getHand, getParam] = parseInput(doHand, g, doParam);

            nBlocks = g.nBlocks();

            obj.order = order;

            % Create the diffOps for each block
            obj.diffOps = cell(1, nBlocks);
            for i = 1:nBlocks
                h = getHand(i);
                p = getParam(i);
                if ~iscell(p)
                    p = {p};
                end
                obj.diffOps{i} = h(g.grids{i}, order, p{:});
            end


            % Build the norm matrix
            H = cell(nBlocks, nBlocks);
            for i = 1:nBlocks
                H{i,i} = obj.diffOps{i}.H;
            end
            obj.H = blockmatrix.toMatrix(H);


            % Build the differentiation matrix
            Ns = zeros(nBlocks,1);
            for i = 1:nBlocks
                Ns(i) = length(obj.diffOps{i}.D);
            end
            obj.blockmatrixDiv = {Ns, Ns};
            D = blockmatrix.zero(obj.blockmatrixDiv);
            for i = 1:nBlocks
                D{i,i} = obj.diffOps{i}.D;
            end

            for i = 1:nBlocks
                for j = 1:nBlocks
                    intf = g.connections{i,j};
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
            obj.grid = g;


            function [getHand, getParam] = parseInput(doHand, g, doParam)
                if ~isa(g, 'multiblock.Grid')
                    error('multiblock:DiffOp:DiffOp:InvalidGrid', 'Requires a multiblock grid.');
                end

                if iscell(doHand) && length(doHand) == g.nBlocks()
                    getHand = @(i)doHand{i};
                elseif isa(doHand, 'function_handle')
                    getHand = @(i)doHand;
                else
                    error('multiblock:DiffOp:DiffOp:InvalidGridDoHand', 'doHand must be a function handle or a cell array of length grid.nBlocks');
                end

                if isempty(doParam)
                    getParam = @(i){};
                    return
                end

                if ~iscell(doParam)
                    getParam = @(i)doParam;
                    return
                end

                % doParam is a non-empty cell-array

                if length(doParam) == g.nBlocks() && all(cellfun(@iscell, doParam))
                    % doParam is a cell-array of cell-arrays
                    getParam = @(i)doParam{i};
                    return
                end

                getParam = @(i)doParam;
            end
        end

        function ops = splitOp(obj, op)
            % Splits a matrix operator into a cell-matrix of matrix operators for
            % each grid.
            ops = sparse2cell(op, obj.NNN);
        end

        % Get a boundary operator specified by opName for the given boundary/BoundaryGroup
        function op = getBoundaryOperator(obj, opName, boundary)
            switch class(boundary)
                case 'cell'
                    localOpName = [opName '_' boundary{2}];
                    blockId = boundary{1};
                    localOp = obj.diffOps{blockId}.(localOpName);

                    div = {obj.blockmatrixDiv{1}, size(localOp,2)};
                    blockOp = blockmatrix.zero(div);
                    blockOp{blockId,1} = localOp;
                    op = blockmatrix.toMatrix(blockOp);
                    return
                case 'multiblock.BoundaryGroup'
                    op = sparse(size(obj.D,1),0);
                    for i = 1:length(boundary)
                        op = [op, obj.getBoundaryOperator(opName, boundary{i})];
                    end
                otherwise
                    error('Unknown boundary indentifier')
            end
        end

        function op = getBoundaryQuadrature(obj, boundary)
            opName = 'H';
            switch class(boundary)
                case 'cell'
                    localOpName = [opName '_' boundary{2}];
                    blockId = boundary{1};
                    op = obj.diffOps{blockId}.(localOpName);

                    return
                case 'multiblock.BoundaryGroup'
                    N = length(boundary);
                    H_bm = cell(N,N);
                    for i = 1:N
                        H_bm{i,i} = obj.getBoundaryQuadrature(boundary{i});
                    end
                    op = blockmatrix.toMatrix(H_bm);
                otherwise
                    error('Unknown boundary indentifier')
            end
        end

        % Creates the closure and penalty matrix for a given boundary condition,
        %    boundary -- the name of the boundary on the form {id,name} where
        %                id is the number of a block and name is the name of a
        %                boundary of that block example: {1,'s'} or {3,'w'}. It
        %                can also be a boundary group
        function [closure, penalty] = boundary_condition(obj, boundary, type)
            switch class(boundary)
                case 'cell'
                    [closure, penalty] = obj.singleBoundaryCondition(boundary, type);
                case 'multiblock.BoundaryGroup'
                    [n,m] = size(obj.D);
                    closure = sparse(n,m);
                    penalty = sparse(n,0);
                    for i = 1:length(boundary)
                        [closurePart, penaltyPart] = obj.boundary_condition(boundary{i}, type);
                        closure = closure + closurePart;
                        penalty = [penalty, penaltyPart];
                    end
                otherwise
                    error('Unknown boundary indentifier')
            end

        end

        function [closure, penalty] = singleBoundaryCondition(obj, boundary, type)
            I = boundary{1};
            name = boundary{2};

            % Get the closure and penaly matrices
            [blockClosure, blockPenalty] = obj.diffOps{I}.boundary_condition(name, type);

            % Expand to matrix for full domain.
            closure = multiblock.local2globalClosure(blockClosure, obj.blockmatrixDiv, I);
            penalty = multiblock.local2globalPenalty(blockPenalty, obj.blockmatrixDiv, I);
        end

        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary)
            error('not implemented')
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
