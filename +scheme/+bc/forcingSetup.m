% Setup the forcing function for the given boundary conditions and data.
% S_sign allows changing the sign of the function to put on different sides in the system of ODEs.
%   default is 1, which the same side as the diffOp.
function S = forcingSetup(diffOp, penalties, bcs, S_sign)
    default_arg('S_sign', 1);

    assertType(bcs, 'cell');
    assertIsMember(S_sign, [1, -1]);

    scheme.bc.verifyFormat(bcs, diffOp);

    [gridData, symbolicData] = parseAndSortData(bcs, penalties, diffOp);

    % Setup penalty function
    O = spzeros(size(diffOp),1);
    function v = S_fun(t)
        v = O;
        for i = 1:length(gridData)
            v = v + gridData{i}.penalty*gridData{i}.func(t);
        end

        for i = 1:length(symbolicData)
            v = v + symbolicData{i}.penalty*symbolicData{i}.func(t, symbolicData{i}.coords{:});
        end

        v = S_sign * v;
    end
    S = @S_fun;
end

% Go through a cell array of boundary condition specifications and return cell arrays
% of structs for grid and symbolic data.
function [gridData, symbolicData] = parseAndSortData(bcs, penalties, diffOp)
    for i = 1:length(bcs)
        [ok, isSymbolic, data] = parseData(bcs{i}, penalties{i}, diffOp.grid)

        if ~ok
            continue % There was no data
        end

        if isSymbolic
            gridData{end+1} = data;
        else
            symbolicData{end+1} = data;
        end
    end
end

function [ok, isSymbolic, dataStruct] = parseData(bc, penalty, grid)
    if ~isfield(bc,'data') || isempty(bc.data)
        isSymbolic = [];
        dataStruct = struct();
        ok = false;
        return
    end
    ok = true;

    nArg = nargin(bc.data);

    if nArg > 1
        % Symbolic data
        isSymbolic = true;
        coord = grid.getBoundary(bc.boundary);
        dataStruct.penalty = penalty;
        dataStruct.func = bc.data;
        dataStruct.coords = num2cell(coord, 1);
    else
        % Grid data
        isSymbolic = false;
        dataStruct.penalty = penalty;
        dataStruct.func = bcs{i}.data;
    end
end
