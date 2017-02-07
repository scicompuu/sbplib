% Creates a matrix of size n,m with the values of val on the diagonals diag.
%   A = stripeMatrix(val,diags,n,m)
function A = stripeMatrix(val,diags,n,m)
    default_arg('m',n);

    D = ones(n,1)*val;
    A = spdiags(D,diags,n,m);
end