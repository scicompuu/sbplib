% Takes a value and returns a string representation of that value.
% If syntaxFlag is true, a string with valid matlab syntax is returned.
function str = toString(value, syntaxFlag)
    default_arg('syntaxFlag',false);

    if syntaxFlag
        error('Not supported yet.')
    end

    str = value2string(value);
end

function str  = value2string(value)
    if isnumeric(value) || ischar(value) || islogical(value)
        str = mat2str(value);
    elseif isstruct(value)
        str = struct2string(value);
    elseif iscell(value)
        str = cell2string(value);
    elseif isa(value,'function_handle')
        str = func2str(value);
    elseif isa(value,'sym')
        str = char(value);
    else
        warning('No string representation for class ''%s''', class(value))
        str = 'NO_STR_REP';
    end
end

function str = cell2string(c)
    if isempty(c)
        str = '{}';
        return
    end

    [n, m] = size(c);

    str = '{';

    for i = 1:n-1
        for j = 1:m-1
            str = [str sprintf('%s, ', value2string(c{i,j}))];
        end
        str = [str sprintf('%s; ', value2string(c{i,end}))];
    end

    for j = 1:m-1
        str = [str sprintf('%s, ', value2string(c{end,j}))];
    end
    str = [str sprintf('%s}', value2string(c{end,end}))];
end

function str = struct2string(s)
    fn = fieldnames(s);

    if length(fn) == 0
        str = '{}';
        return
    end

    str = '{';

    for i = 1:length(fn) - 1
        value = s.(fn{i});
        str = [str sprintf('%s: %s, ',fn{i}, value2string(value))];
    end
    value = s.(fn{end});
    str = [str sprintf('%s: %s}',fn{end}, value2string(value))];
end

