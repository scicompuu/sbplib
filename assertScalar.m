function assertScalar(obj)
    if ~isscalar(obj)
        error('sbplib:assertScalar:notScalar', '"%s" must be scalar, found size "%s"', inputname(1), toString(size(obj)));
    end
end
