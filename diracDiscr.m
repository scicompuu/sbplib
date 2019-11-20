
function d = diracDiscr(g, x_s, m_order, s_order, H)
    % n-dimensional delta function
    % g: cartesian grid
    % x_s: source point coordinate vector, e.g. [x; y] or [x; y; z].
    % m_order: Number of moment conditions
    % s_order: Number of smoothness conditions
    % H: cell array of 1D norm matrices
    assertType(g, 'grid.Cartesian');
    dim = g.d;
    d_1D = cell(dim,1);

    % Allow for non-cell input in 1D
    if dim == 1
        H = {H};
    end
    % Create 1D dirac discr for each coordinate dir.
    for i = 1:dim
        d_1D{i} = diracDiscr1D(x_s(i), g.x{i}, m_order, s_order, H{i});
    end

    d = d_1D{dim};
    for i = dim-1: -1: 1
        % Perform outer product, transpose, and then turn into column vector
        d = (d_1D{i}*d')';
        d = d(:);
    end

end


% Helper function for 1D delta functions
function ret = diracDiscr1D(x_s, x, m_order, s_order, H)

    m = length(x);

    % Return zeros if x0 is outside grid
    if x_s < x(1) || x_s > x(end)
        ret = zeros(size(x));
        return
    else
        tot_order = m_order+s_order; %This is equiv. to the number of equations solved for
        S = [];
        M = [];

        % Get interior grid spacing
        middle = floor(m/2);
        h = x(middle+1) - x(middle); % Use middle point to allow for staggered grids.

        index = sourceIndices(x_s, x, tot_order, h);

        polynomial = (x(index)-x(index(1)))/(x(index(end))-x(index(1)));
        x_0 = (x_s-x(index(1)))/(x(index(end))-x(index(1)));

        quadrature = diag(H);
        quadrature_weights = quadrature(index)/h;

        h_polynomial = polynomial(2)-polynomial(1);
        b = zeros(tot_order,1);

        for i = 1:m_order
            b(i,1) = x_0^(i-1);
        end

        for i = 1:tot_order
            for j = 1:m_order
                M(j,i) = polynomial(i)^(j-1)*h_polynomial*quadrature_weights(i);
            end
        end

        for i = 1:tot_order
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


function I = sourceIndices(x_s, x, tot_order, h)
    % Find the indices that are within range of of the point source location
    I = find(tot_order*h/2 >= abs(x-x_s));

    if length(I) > tot_order
        if length(I) == tot_order + 2
            I = I(2:end-1);
        elseif length(I) == tot_order + 1
            I = I(1:end-1);
        end
    elseif length(I) < tot_order
        if x_s < x(1) + ceil(tot_order/2)*h
            I = 1:tot_order;
        elseif x_s > x(end) - ceil(tot_order/2)*h
            I = length(x)-tot_order+1:length(x);
        else
            if I(end) < length(x)
                I = [I; I(end)+1];
            else
                I = [I(1)-1; I];
            end
        end
    end
end
