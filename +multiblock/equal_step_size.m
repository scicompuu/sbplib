% Choose the number of points on m, M, so that the step size is equal to that for N points on n.
function [M] = equal_step_size(n,N,m)
    M = round((m*(N-1)+n)/n);
end