% Evaluate different function handle for each block in a multiblock.Grid
% Function handles may optionaly take a time argument
% f -- cell array of function handles
%       f{i} = f_i(t,x,y,...)
% t -- optional time point. If not specified, it is assumed that the functions take only spatial arguments.
function gf = evalOn(g, f, t)
    assertType(g, 'multiblock.Grid');
    assertType(f, 'cell');

    default_arg('t', []);

    grids = g.grids;
    nBlocks = length(grids);
    gf = cell(nBlocks, 1);

    if isempty(t)
        for i = 1:nBlocks
            gf{i} = grid.evalOn(grids{i}, f{i});
        end
    else
        dim = nargin(f{1}) - 1;
        for i = 1:nBlocks
            switch dim
                case 1
                    gf{i} = grid.evalOn(grids{i}, @(x)f{i}(t,x));
                case 2
                    gf{i} = grid.evalOn(grids{i}, @(x,y)f{i}(t,x,y));
                case 3
                    gf{i} = grid.evalOn(grids{i}, @(x,y,z)f{i}(t,x,y,z));
            end
        end
    end

    gf = blockmatrix.toMatrix(gf);
end
