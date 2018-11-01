function S = forcingSetup(diffOp, penalties, bcs, S_sign)
    default_arg('S_sign', 1);

    assertType(bcs, 'cell');
    assertIsMember(S_sign, [1, -1]);

    scheme.bc.verifyFormat(bcs, diffOp);

    % % Setup storage arrays
    % closure = spzeros(size(diffOp));
    % gridData = {};
    % symbolicData = {};

    % Loop over bcs and collect data
    for i = 1:length(bcs)
        % [ok, isSym, data] = parseData(bcs{i}, penalties{i}, diffOp.grid)

        % if ~ok
        %     % There was no data
        %     continue
        % end

        % if isSym
        %     gridData{end+1} = data;
        % else
        %     symbolicData{end+1} = data;
        % end
    end


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

function [ok, isSym, dataStruct] = parseData(bc, penalty, grid)
    if ~isfield(bc,'data') || isempty(bc.data)
        ok = false;
        return
    end
    ok = true;

    nArg = nargin(bc.data);

    if nArg > 1
        % Symbolic data
        isSym = true;
        coord = grid.getBoundary(bc.boundary);
        dataStruct.penalty = penalty;
        dataStruct.func = bc.data;
        dataStruct.coords = num2cell(coord, 1);
    else
        % Grid data
        isSym = false;
        dataStruct.penalty = penalty;
        dataStruct.func = bcs{i}.data;
    end
end
