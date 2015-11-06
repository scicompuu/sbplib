function A = spdiag(a,i)
    n = length(a);
    A = spdiags(a,i,n,n);
end