function assertIsMember(v, allowed)
    assert(ismember(v, allowed), 'Expected ''%s'' to be in the set %s', inputname(1), toString(allowed));
end
