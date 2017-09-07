function assertType(obj, type)
    if ~iscell(type)
        if ~isa(obj, type)
            error('sbplib:assertType:wrongType', '"%s" must have type "%s", found "%s"', inputname(1), type, class(obj));
        end
    else
        if ~isAnyOf(obj, type)
            error('sbplib:assertType:wrongType', '"%s" must be one of the types %s, found "%s"', inputname(1), toString(type), class(obj));
        end
    end
end
