% Useful for evaulating forcing functions with different functional expressions for each block
% f: cell array of function handles fi
% f_i = f_i(x1,y,...,t)
% t: time point. If not specified, it is assumed that the functions take only spatial arguments.
function gf = blockEvalOn(g, f, t)
default_arg('t',[]);

grids = g.grids;
nBlocks = length(grids);
gf = cell(nBlocks,1);

if isempty(t)
	for i = 1:nBlocks
		grid.evalOn(grids{i}, f{i} );
	end
else
	dim = nargin(f{1}) - 1;
	for i = 1:nBlocks
		switch dim
		case 1
			gf{i} = grid.evalOn(grids{i}, @(x)f{i}(x,t) );
		case 2
			gf{i} = grid.evalOn(grids{i}, @(x,y)f{i}(x,y,t) );
		case 3
			gf{i} = grid.evalOn(grids{i}, @(x,y,z)f{i}(x,y,z,t) );
		end
	end
end

gf = blockmatrix.toMatrix(gf);