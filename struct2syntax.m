function str = struct2syntax(s)
    warning('Deprecated! Use toString() instead!');
    fn = fieldnames(s);

    if length(fn) == 0
        str = 'struct()';
        return
    end

    str = sprintf('struct(');

    for i = 1:length(fn) - 1
        value = s.(fn{i});
        str = [str sprintf('''%s'', %s, ',fn{i}, valueString(value))];
    end
    value = s.(fn{end});
    str = [str sprintf('''%s'', %s)',fn{end}, valueString(value))];
end

function str  = valueString(value)
    if ischar(value)
        str = ['''' value ''''];
    else
        str = num2str(value);
    end
end