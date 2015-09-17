% Assert that array A has the size s.
function assert_size(A,s)
    errmsg = sprintf('Expected %s to have size %s, got: %s',inputname(1), format_vector(s), format_vector(size(A)));
    assert(all(size(A) == s),errmsg);
end

function str = format_vector(a)
    l = length(a);
    str = sprintf('[%d',a(1));

    for i = 2:l
        str = [str sprintf(', %d',a(i))];
    end

    str = [str ']'];
end