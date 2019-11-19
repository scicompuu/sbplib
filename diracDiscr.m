
function d = diracDiscr(x_s, x, m_order, s_order, H)
    % n-dimensional delta function
    % x_s: source point coordinate vector, e.g. [x, y] or [x, y, z].
    % x: cell array of grid point column vectors for each dimension.
    % m_order: Number of moment conditions
    % s_order: Number of smoothness conditions
    % H: cell array of 1D norm matrices

    dim = length(x_s);
    d_1D = cell(dim,1);

    % If 1D, non-cell input is accepted
    if dim == 1 && ~iscell(x)
        d = diracDiscr1D(x_s, x, m_order, s_order, H);

    else
        for i = 1:dim
            d_1D{i} = diracDiscr1D(x_s(i), x{i}, m_order, s_order, H{i});
        end

        d = d_1D{dim};
        for i = dim-1: -1: 1
            % Perform outer product, transpose, and then turn into column vector
            d = (d_1D{i}*d')';
            d = d(:);
        end
    end

end


% Helper function for 1D delta functions
function ret = diracDiscr1D(x_s , x , m_order, s_order, H)

    m = length(x);

    % Return zeros if x0 is outside grid
    if(x_s < x(1) || x_s > x(end) )

        ret = zeros(size(x));

    else

        fnorm = diag(H);
        tot_order = m_order+s_order; %This is equiv. to the number of equations solved for
        S = [];
        M = [];

        % Get interior grid spacing
        middle = floor(m/2);
        h = x(middle+1) - x(middle);

        % Find the indices that are within range of of the point source location
        ind_delta = find(tot_order*h/2 >= abs(x-x_s));

        % Ensure that ind_delta is not too long
        if length(ind_delta) == (tot_order + 2)
            ind_delta = ind_delta(2:end-1);
        elseif length(ind_delta) == (tot_order + 1)
            ind_delta = ind_delta(1:end-1);
        end

        % Use first tot_order grid points
        if length(ind_delta)<tot_order && x_s < x(1) + ceil(tot_order/2)*h;
            index=1:tot_order;
            polynomial=(x(1:tot_order)-x(1))/(x(tot_order)-x(1));
            x_0=(x_s-x(1))/(x(tot_order)-x(1));
            norm=fnorm(1:tot_order)/h;

        % Use last tot_order grid points
        elseif length(ind_delta)<tot_order && x_s > x(end) - ceil(tot_order/2)*h;
            index = length(x)-tot_order+1:length(x);
            polynomial = (x(end-tot_order+1:end)-x(end-tot_order+1))/(x(end)-x(end-tot_order+1));
            norm = fnorm(end-tot_order+1:end)/h;
            x_0 = (x_s-x(end-tot_order+1))/(x(end)-x(end-tot_order+1));

        % Interior, compensate for round-off errors.
        elseif length(ind_delta) < tot_order
            if ind_delta(end)<m
                ind_delta = [ind_delta; ind_delta(end)+1];
            else
                ind_delta = [ind_delta(1)-1; ind_delta];
            end
            
            index = ind_delta;
            polynomial = (x(ind_delta)-x(ind_delta(1)))/(x(ind_delta(end))-x(ind_delta(1)));
            x_0 = (x_s-x(ind_delta(1)))/(x(ind_delta(end))-x(ind_delta(1)));
            norm = fnorm(ind_delta)/h;

        % Interior
        else
            index = ind_delta;
            polynomial = (x(ind_delta)-x(ind_delta(1)))/(x(ind_delta(end))-x(ind_delta(1)));
            x_0 = (x_s-x(ind_delta(1)))/(x(ind_delta(end))-x(ind_delta(1)));
            norm = fnorm(ind_delta)/h;
        end

        h_polynomial = polynomial(2)-polynomial(1);
        b = zeros(m_order+s_order,1);

        for i = 1:m_order
            b(i,1) = x_0^(i-1);
        end

        for i = 1:(m_order+s_order)
            for j = 1:m_order
                M(j,i) = polynomial(i)^(j-1)*h_polynomial*norm(i);
            end
        end

        for i = 1:(m_order+s_order)
            for j = 1:s_order
                S(j,i) = (-1)^(i-1)*polynomial(i)^(j-1);
            end
        end

        A = [M;S];

        d = A\b;
        ret = x*0;
        ret(index) = d/h*h_polynomial;
    end

end






