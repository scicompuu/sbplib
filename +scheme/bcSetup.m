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
function [closure, S] = bcSetup(diffOp, bc, S_sign)
    default_arg('S_sign', 1);
    assertType(bc, 'cell');
    assert(S_sign == 1 || S_sign == -1, 'S_sign must be either 1 or -1');


    % Setup storage arrays
    closure = spzeros(size(diffOp));
    gridDataPenalties = {};
    gridDataFunctions = {};
    symbolicDataPenalties = {};
    symbolicDataFunctions = {};
    symbolicDataCoords = {};

    % Collect closures, penalties and data
    for i = 1:length(bc)
        assertType(bc{i}, 'struct');
        [localClosure, penalty] = diffOp.boundary_condition(bc{i}.boundary, bc{i}.type);
        closure = closure + localClosure;

        if ~isfield(bc{i},'data') || isempty(bc{i}.data)
            % Skip to next loop if there is no data
            continue
        end
        assertType(bc{i}.data, 'function_handle');

        % Find dimension
        dim = size(diffOp.grid.getBoundary(bc{i}.boundary), 2);

        if nargin(bc{i}.data) == 1
            % Grid data
            boundarySize = [size(diffOp.grid.getBoundary(bc{i}.boundary),1),1];
            assert_size(bc{i}.data(0), boundarySize); % Eval for t = 0 and make sure the function returns a grid vector of the correct size.
            gridDataPenalties{end+1} = penalty;
            gridDataFunctions{end+1} = bc{i}.data;
        elseif nargin(bc{i}.data) == 1+dim
            % Symbolic data
            coord = diffOp.grid.getBoundary(bc{i}.boundary);
            symbolicDataPenalties{end+1} = penalty;
            symbolicDataFunctions{end+1} = bc{i}.data;
            symbolicDataCoords{end+1} = num2cell(coord ,1);
        else
            error('sbplib:scheme:bcSetup:DataWrongNumberOfArguments', 'bc{%d}.data has the wrong number of input arguments. Must be either only time or time and space.', i);
        end
    end

    % Setup penalty function
    O = spzeros(size(diffOp),1);
    function v = S_fun(t)
        v = O;
        for i = 1:length(gridDataFunctions)
            v = v + gridDataPenalties{i}*gridDataFunctions{i}(t);
        end

        for i = 1:length(symbolicDataFunctions)
            v = v + symbolicDataPenalties{i}*symbolicDataFunctions{i}(t, symbolicDataCoords{i}{:});
        end

        v = S_sign * v;
    end
    S = @S_fun;
end
