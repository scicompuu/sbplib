function A = spdiag(a,i)
    default_arg('i',0);

    if isrow(a)
        a = a';
    end

    n = length(a)-abs(i);
    A = spdiags(a,i,n,n);
end