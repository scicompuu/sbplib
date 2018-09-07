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
