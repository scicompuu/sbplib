function A = spdiag(a,i)
    n = length(a)-abs(i);
    A = spdiags(a,i,n,n);
end