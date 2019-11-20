function d = diracDiscrCurve(x_s, g, m_order, s_order, order, opSet)
    % 2-dimensional delta function for single-block curvilinear grid
    % x_s:      source point coordinate vector, e.g. [x; y] or [x; y; z].
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

    % Compute Jacobian
    J = jacobian(g, opSet, order);

    % Find approximate logical coordinates of point source
    X = g.points();
    U = g.logic.points();
    U_interp = scatteredInterpolant(X, U(:,1));
    V_interp = scatteredInterpolant(X, U(:,2));
    uS = U_interp(x_s);
    vS = V_interp(x_s);

    % Get quadrature matrices for moment conditions
    m = g.size();
    ops_u = opSet{1}(m(1), {0, 1}, order);
    ops_v = opSet{2}(m(2), {0, 1}, order);
    H_u =  ops_u.H;
    H_v =  ops_v.H;

    % Get delta function for logical grid and scale by Jacobian
    d = (1./J).*diracDiscr(g, [uS; vS], m_order, s_order, {H_u, H_v});

end

function J = jacobian(g, opSet, order)
    m = g.size();
    m_u = m(1);
    m_v = m(2);
    ops_u = opSet{1}(m_u, {0, 1}, order);
    ops_v = opSet{2}(m_v, {0, 1}, order);
    I_u = speye(m_u);
    I_v = speye(m_v);

    D1_u = ops_u.D1;
    D1_v = ops_v.D1;

    Du = kr(D1_u,I_v);
    Dv = kr(I_u,D1_v);

    coords = g.points();
    x = coords(:,1);
    y = coords(:,2);

    x_u = Du*x;
    x_v = Dv*x;
    y_u = Du*y;
    y_v = Dv*y;

    J = x_u.*y_v - x_v.*y_u;
end