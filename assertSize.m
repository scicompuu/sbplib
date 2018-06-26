% Assert that array A has the size s.
function assertSize(A,s)
    errmsg = sprintf('Expected %s to have size %s, got: %s',inputname(1), toString(s), toString(size(A)));
    assert(all(size(A) == s),errmsg);
end
