function str = struct2string(s)
    warning('Deprecated! Use toString() instead!');
    fn = fieldnames(s);

    if length(fn) == 0
        str = '{}';
        return
    end

    str = sprintf('{');

    for i = 1:length(fn) - 1
        value = s.(fn{i});
        str = [str sprintf('%s: %s, ',fn{i}, valueString(value))];
    end
    value = s.(fn{end});
    str = [str sprintf('%s: %s}',fn{end}, valueString(value))];
end

function str  = value2string(value)
    if isnumeric(value) || ischar(value)
        str = mat2str(value);
    elseif isstruct(value)
        str = struct2string(value);
    else
        str = 'NO_STR_REP';
    end
end