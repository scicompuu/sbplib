% Returns true of obj is any of he types in the cell array types
%    b = isAnyOf(obj, types)
function b = isAnyOf(obj, types)
    for i = 1:length(types)
        if isa(obj, types{i});
            b = true;
            return
        end
    end
    b = false;
end