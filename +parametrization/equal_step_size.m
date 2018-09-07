% Calculates M so that the stepsize m/M is as close to n/M as possible
function M = equal_step_size(n,N,m)
    M = round((m*(N-1)+n)/n);
end