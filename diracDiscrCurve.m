function d = diracDiscrCurve(x_s, g, m_order, s_order, order, opSet)
    % 2-dimensional delta function for single-block curvilinear grid
    % x_s:      source point coordinate vector, e.g. [x, y] or [x, y, z].
    % g:        single-block grid containing the source
    % m_order:  Number of moment conditions
    % s_order:  Number of smoothness conditions
    % order:    Order of SBP derivative approximations
    % opSet:    Cell array of function handle to opSet generator

    default_arg('order', m_order);
    default_arg('opSet', {@sbp.D2Variable, @sbp.D2Variable});

    dim = length(x_s);
    assert(dim == 2, 'diracDiscrCurve: Only implemented for 2d.');
    assert(isa(g, 'grid.Curvilinear'));

    m = g.size();
    m_u = m(1);
    m_v = m(2);
    ops_u = opSet{1}(m_u, {0, 1}, order);
    ops_v = opSet{2}(m_v, {0, 1}, order);
    I_u = speye(m_u);
    I_v = speye(m_v);

    D1_u = ops_u.D1;
    H_u =  ops_u.H;

    D1_v = ops_v.D1;
    H_v =  ops_v.H;

    Du = kr(D1_u,I_v);
    Dv = kr(I_u,D1_v);

    u = ops_u.x;
    v = ops_v.x;

    % Compute Jacobian
    coords = g.points();
    x = coords(:,1);
    y = coords(:,2);

    x_u = Du*x;
    x_v = Dv*x;
    y_u = Du*y;
    y_v = Dv*y;

    J = x_u.*y_v - x_v.*y_u;

    % Find approximate logical coordinates of point source
    [U, V] = meshgrid(u, v);
    U_interp = scatteredInterpolant(coords, U(:));
    V_interp = scatteredInterpolant(coords, V(:));
    uS = U_interp(x_s);
    vS = V_interp(x_s);

    d = (1./J).*diracDiscr([uS, vS], {u, v}, m_order, s_order, {H_u, H_v});

end