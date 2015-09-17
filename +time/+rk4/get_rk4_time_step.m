% Calculate the size of the largest time step given the largest evalue for a operator with pure imaginary e.values.
function k = get_rk4_time_step(lambda,l_type)
    default_arg('l_type','complex')

    rad = abs(lambda);
    if strcmp(l_type,'real')
        % Real eigenvalue
        % kl > -2.7852
        k = 2.7852/rad;

    elseif strcmp(l_type,'imag')
        % Imaginary eigenvalue
        % |kl| < 2.8284
        k = 2.8284/rad;
    elseif strcmp(l_type,'complex')
        % |kl| < 2.5
        k = 2.5/rad;
    else
        error('l_type must be one of ''real'',''imag'' or ''complex''.')
    end
end