function A = spdiag(a)
    n = length(a);
    A = spdiags(a,0,n,n);
end