% Creates a periodic discretization matrix of size n x n
%  with the values of val on the diagonals diag.
%   A = stripeMatrix(val,diags,n)
function A = stripeMatrixPeriodic(val,diags,n)

    D = ones(n,1)*val;
    A = spdiagsPeriodic(D,diags);
end