function varargout = dealStruct(s, fields)
    default_arg('fields', []);

    if isempty(fields)
        out = dealFields(s, fieldnames(s));
        varargout = out(1:nargout);
    else
        assert(nargout == length(fields), 'Number of output arguements must match the number of fieldnames provided');
        varargout = dealFields(s, fields);
    end
end

function out = dealFields(s, fields)
    out = cell(1, length(fields));
    for i = 1:length(fields)
        out{i} = s.(fields{i});
    end
end
