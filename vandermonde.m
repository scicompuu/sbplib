% Create vandermonde matrix for points x and polynomials of order p
% x and p are vectors
% v is a length(x) by length(p) matrix
function V = vandermonde(x, p)
    V = sym(zeros(length(x), length(p))); % Is there a way to make this work for both double and sym

    for i = 1:length(p)
        V(:, i) = mononomial(x,p(i));
    end
end
