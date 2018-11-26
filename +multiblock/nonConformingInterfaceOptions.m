% g: 			multiblock grid
% order:		cell array of the orders of accuracy used in the different blocks,
%				or a scalar for the same order everywhere.
% interpOpSet:  string, e.g 'MC' or 'AWW' The same interpOpSet is used everywhere.
%
% Returns an InterfaceOptions object that can be used as the 'interfaceOptions' argument to multiblock.DiffOp
function options = nonConformingInterfaceOptions(g, orders, interpOpSet)
default_arg(interpOpSet, 'AWW');
nBlocks = g.nBlocks;
conn = g.connections;

% If order is a scalar, the same order is used in all blocks
if ~iscell(orders)
	o = orders;
	orders = cell(1,nBlocks);
	for i = 1:nBlocks
		orders{i} = o;
	end
end

interpOpts = cell(nBlocks, nBlocks);
for i = 1:nBlocks
	for j = 1:nBlocks
		intf = conn{i,j};
        if isempty(intf)
            continue
        end

        mi = length( g.getBoundary({i, intf{1}}) ) - 1;
    	mj = length( g.getBoundary({j, intf{2}}) ) - 1;

		if mi == mj
			% Matching grids, no interpolation required (presumably)
			continue;
		elseif mi/mj == 2
			% Block i is finer

			switch interpOpSet
			case 'MC'
				interpOpSet = sbp.InterpMC(mj+1, mi+1, orders{j}, orders{i});
				I_i2j_good = interpOpSet.IF2C;
                I_i2j_bad = interpOpSet.IF2C;
                I_j2i_good = interpOpSet.IC2F;
                I_j2i_bad = interpOpSet.IC2F;

            case 'AWW'
            	interpOpSetF2C = sbp.InterpAWW(mj+1, mi+1, orders{j}, orders{i}, 'F2C');
            	interpOpSetC2F = sbp.InterpAWW(mj+1, mi+1, orders{j}, orders{i}, 'C2F');
				I_i2j_good = interpOpSetF2C.IF2C;
                I_i2j_bad = interpOpSetC2F.IF2C;
                I_j2i_good = interpOpSetC2F.IC2F;
                I_j2i_bad = interpOpSetF2C.IC2F;
            end

		elseif mj/mi == 2
			% Block j is finer

			switch interpOpSet
			case 'MC'
				interpOpSet = sbp.InterpMC(mi+1, mj+1, orders{i}, orders{j});
				I_i2j_good = interpOpSet.IC2F;
                I_i2j_bad = interpOpSet.IC2F;
                I_j2i_good = interpOpSet.IF2C;
                I_j2i_bad = interpOpSet.IF2C;

            case 'AWW'
            	interpOpSetF2C = sbp.InterpAWW(mi+1, mj+1, orders{i}, orders{j}, 'F2C');
            	interpOpSetC2F = sbp.InterpAWW(mi+1, mj+1, orders{i}, orders{j}, 'C2F');
				I_i2j_good = interpOpSetC2F.IC2F;
                I_i2j_bad = interpOpSetF2C.IC2F;
                I_j2i_good = interpOpSetF2C.IF2C;
                I_j2i_bad = interpOpSetC2F.IF2C;
            end
		else
			error(sprintf('Interpolation operators for grid ratio %f have not yet been constructed', mi/mj));
		end

		interpOpts{i,j} = cell(2,1);
		interpOpts{i,j}{1}.I_local2neighbor.good = I_i2j_good;
		interpOpts{i,j}{1}.I_local2neighbor.bad = I_i2j_bad;
		interpOpts{i,j}{1}.I_neighbor2local.good = I_j2i_good;
		interpOpts{i,j}{1}.I_neighbor2local.bad = I_j2i_bad;

		interpOpts{i,j}{2}.I_local2neighbor.good = I_j2i_good;
		interpOpts{i,j}{2}.I_local2neighbor.bad = I_j2i_bad;
		interpOpts{i,j}{2}.I_neighbor2local.good = I_i2j_good;
		interpOpts{i,j}{2}.I_neighbor2local.bad = I_i2j_bad;

		options = multiblock.InterfaceOptions(g, interpOpts);

	end
end

