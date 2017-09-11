% Tests if consecutive elements of vector v are euidistant
function b = isEquidistant(v)
    if length(v) < 2
        error('sbplib:isEquidistant:inputTooShort', 'Input vector is too short');
    end

    tol = 1e-8;

    d = v(2:end) - v(1:end-1);
    err = abs(d - d(1));

    relErr = err./abs(d);

    I_zero = find(d < tol);

    relErr(I_zero) = err(I_zero);

    b = all(relErr < tol);
end
