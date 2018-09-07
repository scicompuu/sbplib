% function [closure, S] = bcSetup(diffOp, bc)
% Takes a diffOp and a cell array of boundary condition definitions.
% Each bc is a struct with the fields
%  * type     -- Type of boundary condition
%  * boundary -- Boundary identifier
%  * data     -- A function_handle for a function which provides boundary data.(see below)
% Also takes S_sign which modifies the sign of S, [-1,1]
% Returns a closure matrix and a forcing function S.
%
% The boundary data function can either be a function of time or a function of time and space coordinates.
% In the case where it only depends on time it should return the data as grid function for the boundary.
% In the case where it also takes space coordinates the number of space coordinates should match the number of dimensions of the problem domain.
% For example in the 2D case: f(t,x,y).
function [closure, S] = bcSetup(diffOp, bcs, S_sign)
    default_arg('S_sign', 1);
    assertType(bcs, 'cell');
    assert(S_sign == 1 || S_sign == -1, 'S_sign must be either 1 or -1');

    verifyBcFormat(bcs, diffOp);

    % Setup storage arrays
    closure = spzeros(size(diffOp));
    gridData = {};
    symbolicData = {};

    % Collect closures, penalties and data
    for i = 1:length(bcs)
        [localClosure, penalty] = diffOp.boundary_condition(bcs{i}.boundary, bcs{i}.type);
        closure = closure + localClosure;

        [ok, isSym, data] = parseData(bcs{i}, penalty, diffOp.grid);

        if ~ok
            % There was no data
            continue
        end

        if isSym
            symbolicData{end+1} = data;
        else
            gridData{end+1} = data;
        end
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

function verifyBcFormat(bcs, diffOp)
    for i = 1:length(bcs)
        assertType(bcs{i}, 'struct');
        assertStructFields(bcs{i}, {'type', 'boundary'});

        if ~isfield(bcs{i}, 'data') || isempty(bcs{i}.data)
            continue
        end

        if ~isa(bcs{i}.data, 'function_handle')
            error('bcs{%d}.data should be a function of time or a function of time and space',i);
        end

        b = diffOp.grid.getBoundary(bcs{i}.boundary);

        dim = size(b,2);

        if nargin(bcs{i}.data) == 1
            % Grid data (only function of time)
            assertSize(bcs{i}.data(0), 1, size(b));
        elseif nargin(bcs{i}.data) ~= 1+dim
           error('sbplib:scheme:bcSetup:DataWrongNumberOfArguments', 'bcs{%d}.data has the wrong number of input arguments. Must be either only time or time and space.', i);
        end
    end
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
