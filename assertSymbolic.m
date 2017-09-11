function assertSymbolic(s)
    assert(logical(simplify(s)));
end
