function default_struct(s, val)
    if evalin('caller',sprintf('~exist(''%s'',''var'')',s))
        given = [];
    else
        given = evalin('caller', s);
    end

    final = copyWithDefault(given, val);
    assignin('caller', s, final);
end
