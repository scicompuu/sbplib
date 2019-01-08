function assertLength(A,l)
    assert(isvector(A), sprintf('Expected ''%s'' to be a vector, got matrix of size %s',inputname(1), toString(size(A))));
    assert(length(A) == l, sprintf('Expected ''%s'' to have length %d, got %d', inputname(1), l, length(A)));
end
