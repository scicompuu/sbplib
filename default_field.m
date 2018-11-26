function default_field(s, f, val)
    if isfield(s,f)
    	field = getfield(s, f);
    	if ~isempty(field)
        	return
        end
    end
    s.(f) = val;
    assignin('caller', inputname(1),s);
end