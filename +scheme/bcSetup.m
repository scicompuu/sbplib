% function [closure, S] = bcSetup(diffOp, bc)
% Takes a diffOp and a cell array of boundary condition definitions.
% Each bc is a struct with the fields
%  * type     -- Type of boundary condition
%  * boundary -- Boundary identifier
%  * data     -- A function_handle with time and space coordinates as a parameters, for example f(t,x,y) for a 2D problem
% Also takes S_sign which modifies the sign of S, [-1,1]
% Returns a closure matrix and a forcing function S
function [closure, S] = bcSetup(diffOp, bc, S_sign)
    default_arg('S_sign', 1);
    assertType(bc, 'cell');
    assert(S_sign == 1 || S_sign == -1, 'S_sign must be either 1 or -1');


    closure = spzeros(size(diffOp));
    penalties = {};
    dataFunctions = {};
    dataParams = {};

    for i = 1:length(bc)
        assertType(bc{i}, 'struct');
        [localClosure, penalty] = diffOp.boundary_condition(bc{i}.boundary, bc{i}.type);
        closure = closure + localClosure;

        if ~isfield(bc{i},'data') || isempty(bc{i}.data)
            continue
        end
        assertType(bc{i}.data, 'function_handle');

        coord = diffOp.grid.getBoundary(bc{i}.boundary);
        assertNumberOfArguments(bc{i}.data, 1+size(coord,2));

        penalties{end+1} = penalty;
        dataFunctions{end+1} = bc{i}.data;
        dataParams{end+1} = num2cell(coord ,1);
    end

    O = spzeros(size(diffOp),1);
    function v = S_fun(t)
        v = O;
        for i = 1:length(dataFunctions)
            v = v + penalties{i}*dataFunctions{i}(t, dataParams{i}{:});
        end

        v = S_sign * v;
    end
    S = @S_fun;
end
