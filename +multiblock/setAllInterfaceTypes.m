% Create interface configuration with a single type for all interfaces
%    g -- multiblock grid
% type -- type for all interfaces
function intfTypes = setAllInterfaceTypes(g, type)
	intfTypes = cell(g.nBlocks(), g.nBlocks());
	for i = 1:g.nBlocks()
		for j = 1:g.nBlocks()
			if ~isempty(g.connections{i,j})
				intfTypes{i,j} = type;
			end
		end
	end
end
