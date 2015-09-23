% Returns true if the the fields of struct a exists in A and have the same values
function b = structIsSubset(a,A)
    fn = fieldnames(a);

    b = true; % if a has no filds
    for j = 1:length(fn)
        fname = fn{j};
        value = a.(fname);
        if isfield(A,fname) && a.(fname) == A.(fname)
            b = true;
            continue;
        else
            b = false;
            break;
        end
    end
end