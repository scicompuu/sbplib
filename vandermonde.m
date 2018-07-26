% Create vandermonde matrix for points x and polynomials of order p
% x is a list of N points of size [N,dim],
% p is a list of polynomial orders of size [M, dim].
% the given mononomials are evaluated and the NxM matrix V is returned.
function V = vandermonde(x, p)
    assert(size(x,2) == size(p,2), 'x and p must have the same number of columns')
    n = size(x,1);
    m = size(p,1);

    for i = 1:m
        V(:,i) = mononomial(x, p(i,:));
    end

    assertSize(V,[n,m]);
end
