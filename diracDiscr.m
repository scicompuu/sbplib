
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
function ret = diracDiscr1D(x_0in , x , m_order, s_order, H)

m = length(x);

% Return zeros if x0 is outside grid
if(x_0in < x(1) || x_0in > x(end) )

    ret = zeros(size(x));

else

    fnorm = diag(H);
    eta = abs(x-x_0in);
    tot = m_order+s_order;
    S = [];
    M = [];

    % Get interior grid spacing
    middle = floor(m/2);
    h = x(middle+1) - x(middle);

    poss = find(tot*h/2 >= eta);

    % Ensure that poss is not too long
    if length(poss) == (tot + 2)
        poss = poss(2:end-1);
    elseif length(poss) == (tot + 1)
        poss = poss(1:end-1);
    end

    % Use first tot grid points
    if length(poss)<tot && x_0in < x(1) + ceil(tot/2)*h;
        index=1:tot;
        pol=(x(1:tot)-x(1))/(x(tot)-x(1));
        x_0=(x_0in-x(1))/(x(tot)-x(1));
        norm=fnorm(1:tot)/h;

    % Use last tot grid points
    elseif length(poss)<tot && x_0in > x(end) - ceil(tot/2)*h;
        index = length(x)-tot+1:length(x);
        pol = (x(end-tot+1:end)-x(end-tot+1))/(x(end)-x(end-tot+1));
        norm = fnorm(end-tot+1:end)/h;
        x_0 = (x_0in-x(end-tot+1))/(x(end)-x(end-tot+1));

    % Interior, compensate for round-off errors.
    elseif length(poss) < tot
        if poss(end)<m
            poss = [poss; poss(end)+1];
        else
            poss = [poss(1)-1; poss];
        end
        pol = (x(poss)-x(poss(1)))/(x(poss(end))-x(poss(1)));
        x_0 = (x_0in-x(poss(1)))/(x(poss(end))-x(poss(1)));
        norm = fnorm(poss)/h;
        index = poss;

    % Interior
    else
        pol = (x(poss)-x(poss(1)))/(x(poss(end))-x(poss(1)));
        x_0 = (x_0in-x(poss(1)))/(x(poss(end))-x(poss(1)));
        norm = fnorm(poss)/h;
        index = poss;
    end

    h_pol = pol(2)-pol(1);
    b = zeros(m_order+s_order,1);

    for i = 1:m_order
        b(i,1) = x_0^(i-1);
    end

    for i = 1:(m_order+s_order)
        for j = 1:m_order
            M(j,i) = pol(i)^(j-1)*h_pol*norm(i);
        end
    end

    for i = 1:(m_order+s_order)
        for j = 1:s_order
            S(j,i) = (-1)^(i-1)*pol(i)^(j-1);
        end
    end

    A = [M;S];

    d = A\b;
    ret = x*0;
    ret(index) = d/h*h_pol;
end

end






