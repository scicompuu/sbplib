% Find the equation for the stencil for d^k/dx^k
function [A,b] = stencilEquation(k, offsets, order)
    q = sym('q', [1, length(offsets)]);

    p = 0:(order-1+k);

    v     = vandermonde(offsets, p);
    vdiff = vandermonde(      0, p-k);

    eq = q*v == vdiff;

    [A,b] = equationsToMatrix(eq, q);
end
