% Errors with a more or less detailed error message if there is a problem with the bc specification
function verifyBcFormat(bcs, diffOp)
    assertType(bcs, 'cell');
    for i = 1:length(bcs)
        assertType(bcs{i}, 'struct');
        assertStructFields(bcs{i}, {'type', 'boundary'});

        if ~isfield(bcs{i}, 'data') || isempty(bcs{i}.data)
            continue
        end

        if ~isa(bcs{i}.data, 'function_handle')
            error('bcs{%d}.data should be a function of time or a function of time and space',i);
        end

        % Find dimension of boundary
        b = diffOp.grid.getBoundary(bcs{i}.boundary);
        dim = size(b,2);

        % Assert that the data function has a valid number of input arguments
        if ~(nargin(bcs{i}.data) == 1 || nargin(bcs{i}.data) == 1 + dim)
            error('sbplib:scheme:bcSetup:DataWrongNumberOfArguments', 'bcs{%d}.data has the wrong number of input arguments. Must be either only time or time and space.', i);
        end

        if nargin(bcs{i}.data) == 1
            % Grid data (only function of time)
            % Assert that the data has the correct dimension
            assertSize(bcs{i}.data(0), 1, size(b));
        end
    end
end
