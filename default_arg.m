function default_arg(s, val)
    if evalin('caller',sprintf('~exist(''%s'',''var'') || isempty(%s)',s,s))
        assignin('caller',s,val)
    end
end
