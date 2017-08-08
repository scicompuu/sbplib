function assertType(obj, type)
    if ~isa(obj, type)
        error('sbplib:assertType:wrongType', '"%s" must have type "%s", found "%s"', inputname(1), type, class(obj));
    end
end
