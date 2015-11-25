function default_field(s, f, val)
    if isfield(s,f)
        return
    end
    s.(f) = val;
    assignin('caller', inputname(1),s);
end