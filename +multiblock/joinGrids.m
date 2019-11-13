% Connects several multiblock grids into one grid
% gs 		   -- 	Cell array of multiblock grids
% gConnections -- 	Upper-triangular cell matrix
% 					gConnections{i,j} specifies all connections between grid i and grid j
%					Example:
%					gConnections{i,j} = { {{1,'e'},{2,'w'}}, {{5,'s'},{2,'n'}},... };
% names 	   --	(Optional) cell array of strings, used for boundary groups in new grid
%					default: names = {'g1', 'g2', ..., 'gN'};
% 					Boundary groups from grid i are contained in g.boundaryGroups.(names{i}), etc.
function g = joinGrids(gs, gConnections, names)

	nGrids = numel(gs);

	% Default names are {'g1', 'g2', ... 'gN'}.
	defaultNames = cell(nGrids, 1);
	for i = 1:nGrids
		defaultNames{i} = sprintf('g%d',i);
	end
	default_arg('names', defaultNames);

	nBlocks = 0;
	for i = 1:nGrids
		nBlocks = nBlocks + gs{i}.nBlocks();
	end

	% Create vector of cumulative sum of number of blocks per grid
	startIndex = zeros(1, nGrids);
	for i = 2:nGrids
		startIndex(i) = startIndex(i-1) + gs{i-1}.nBlocks();
	end

	% Create cell array of all grids
	grids = cell(nBlocks, 1);
	for i = 1:nGrids
		for j = 1:gs{i}.nBlocks();
			grids{startIndex(i)+j} = gs{i}.grids{j};
		end
	end

	% Create cell matrix of connections
	connections = cell(nBlocks, nBlocks);

	% Connections within grids
	for i = 1:nGrids
		for j = 1:gs{i}.nBlocks()
			for k = 1:gs{i}.nBlocks()
				connections{startIndex(i)+j,startIndex(i)+k} = gs{i}.connections{j,k};
			end
		end
	end

	% Connections between grids
	for i = 1:nGrids
		for j = 1:nGrids
			for k = 1:numel(gConnections{i,j})
				b1 = gConnections{i,j}{k}{1};
				id1 = b1{1};
				str1 = b1{2};

				b2 = gConnections{i,j}{k}{2};
				id2 = b2{1};
				str2 = b2{2};

				connections{startIndex(i)+id1, startIndex(j)+id2} = {str1, str2};
			end
		end
	end

	% Boundary groups
	boundaryGroups = struct;
	for i = 1:nGrids
		bgs = gs{i}.boundaryGroups;
		bgNames = fieldnames(bgs);
		for j = 1:numel(bgNames)
			bg = bgs.(bgNames{j});

			% Shift block id:s in boundary groups
			for k = 1:length(bg)
				bg{k}{1} = bg{k}{1} + startIndex(i);
			end

			bgs.(bgNames{j}) = bg;
		end
		boundaryGroups.(names{i}) = bgs;
	end

	% Create grid object
	g = multiblock.Grid(grids, connections, boundaryGroups);

end