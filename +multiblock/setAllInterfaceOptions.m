% g: multiblock grid
% opts: struct of options
% Returns g.nBlocks x g.nBlocks cell array with all elements set to opts.
function optsCell = setAllInterfaceOptions(g, opts)

	optsCell = cell(g.nBlocks(), g.nBlocks());
	for i = 1:g.nBlocks()^2
		optsCell{i} = opts;
	end

end