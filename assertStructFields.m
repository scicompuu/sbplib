% Assert that the struct s has the all the field names in the cell array fns.
function assertStructFields(s, fns)
    assertType(s, 'struct');
    assertType(fns, 'cell');

    ok = ismember(fns, fieldnames(s));
    if ~all(ok)
        str1 = sprintf("'%s' must have the fields %s\n", inputname(1), toString(fns));
        str2 = sprintf("The following fields are missing: %s", toString(fns(~ok)));
        error(str1 + str2);
    end
end
