function varargout = dealStruct(s, fields)
    default_arg('fields', fieldnames(s));

    assert(nargout == length(fields), 'Number of output arguements must match the number of fields');

    for i = 1:length(fields)
        varargout{i} = s.(fields{i});
    end
end
