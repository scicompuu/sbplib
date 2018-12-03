% g: multiblock grid
% opts: struct of options
% Returns g.nBlocks x g.nBlocks cell array with all interface opts set to opts.
function optsCell = setAllInterfaceTypes(g, opts)

	optsCell = cell(g.nBlocks(), g.nBlocks());
	for i = 1:g.nBlocks()
		for j = 1:g.nBlocks
			if ~isempty(g.connections{i,j})
				optsCell{i,j} = opts;
			end
		end
	end

end