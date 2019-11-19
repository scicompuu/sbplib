function tests = diracDiscrTest()
	    tests = functiontests(localfunctions);
end

%TODO: Test discretizing with smoothness conditions.
%      Only discretization with moment conditions currently tested.

function testLeftRandom(testCase)

    orders = [2, 4, 6];
    mom_conds = orders;
    rng(0xDABBAD00) % Set seed. Jabba-dabba-doooo!

    for o = 1:length(orders)
        order = orders(o);
        mom_cond = mom_conds(o);
        [xl, ~, h, x, H, fs] = setup1D(order, mom_cond);

        % Test random points near left boundary
        x0s = xl + 2*h*rand(1,10);

        for j = 1:length(fs)
                f = fs{j};
                fx = f(x);
            for i = 1:length(x0s)
                x0 = x0s(i);
                delta = diracDiscr(x0, x, mom_cond, 0, H);
                integral = delta'*H*fx;
                err = abs(integral - f(x0));
                testCase.verifyLessThan(err, 1e-12);
            end
        end
    end
end

function testRightRandom(testCase)

    orders = [2, 4, 6];
    mom_conds = orders;
    rng(0xDABBAD00) % Set seed. Jabba-dabba-doooo!

    for o = 1:length(orders)
        order = orders(o);
        mom_cond = mom_conds(o);
        [~, xr, h, x, H, fs] = setup1D(order, mom_cond);

        % Test random points near right boundary
        x0s = xr - 2*h*rand(1,10);

        for j = 1:length(fs)
                f = fs{j};
                fx = f(x);
            for i = 1:length(x0s)
                x0 = x0s(i);
                delta = diracDiscr(x0, x, mom_cond, 0, H);
                integral = delta'*H*fx;
                err = abs(integral - f(x0));
                testCase.verifyLessThan(err, 1e-12);
            end
        end
    end
end

function testInteriorRandom(testCase)

    orders = [2, 4, 6];
    mom_conds = orders;
    rng(0xDABBAD00) % Set seed. Jabba-dabba-doooo!

    for o = 1:length(orders)
        order = orders(o);
        mom_cond = mom_conds(o);
        [xl, xr, h, x, H, fs] = setup1D(order, mom_cond);

        % Test random points in interior
        x0s = (xl+2*h) + (xr-xl-4*h)*rand(1,20);

        for j = 1:length(fs)
                f = fs{j};
                fx = f(x);
            for i = 1:length(x0s)
                x0 = x0s(i);
                delta = diracDiscr(x0, x, mom_cond, 0, H);
                integral = delta'*H*fx;
                err = abs(integral - f(x0));
                testCase.verifyLessThan(err, 1e-12);
            end
        end
    end
end

% x0 outside grid should yield 0 integral!
function testX0OutsideGrid(testCase)

    orders = [2, 4, 6];
    mom_conds = orders;

    for o = 1:length(orders)
        order = orders(o);
        mom_cond = mom_conds(o);
        [xl, xr, h, x, H, fs] = setup1D(order, mom_cond);

        % Test points outisde grid
        x0s = [xl-1.1*h, xr+1.1*h];

        for j = 1:length(fs)
                f = fs{j};
                fx = f(x);
            for i = 1:length(x0s)
                x0 = x0s(i);
                delta = diracDiscr(x0, x, mom_cond, 0, H);
                integral = delta'*H*fx;
                err = abs(integral - 0);
                testCase.verifyLessThan(err, 1e-12);
            end
        end
    end
end

function testAllGP(testCase)

    orders = [2, 4, 6];
    mom_conds = orders;

    for o = 1:length(orders)
        order = orders(o);
        mom_cond = mom_conds(o);
        [~, ~, ~, x, H, fs] = setup1D(order, mom_cond);

        % Test all grid points
        x0s = x;

        for j = 1:length(fs)
                f = fs{j};
                fx = f(x);
            for i = 1:length(x0s)
                x0 = x0s(i);
                delta = diracDiscr(x0, x, mom_cond, 0, H);
                integral = delta'*H*fx;
                err = abs(integral - f(x0));
                testCase.verifyLessThan(err, 1e-12);
            end
        end
    end
end

function testHalfGP(testCase)

    orders = [2, 4, 6];
    mom_conds = orders;

    for o = 1:length(orders)
        order = orders(o);
        mom_cond = mom_conds(o);
        [~, ~, ~, x, H, fs] = setup1D(order, mom_cond);

        % Test halfway between all grid points
        x0s = 1/2*( x(2:end)+x(1:end-1) );

        for j = 1:length(fs)
                f = fs{j};
                fx = f(x);
            for i = 1:length(x0s)
                x0 = x0s(i);
                delta = diracDiscr(x0, x, mom_cond, 0, H);
                integral = delta'*H*fx;
                err = abs(integral - f(x0));
                testCase.verifyLessThan(err, 1e-12);
            end
        end
    end
end

%=============== 2D tests ==============================
function testAllGP2D(testCase)

    orders = [2, 4, 6];
    mom_conds = orders;

    for o = 1:length(orders)
        order = orders(o);
        mom_cond = mom_conds(o);
        [~, ~, x, X, ~, H, fs] = setup2D(order, mom_cond);
        H_global = kron(H{1}, H{2});

        % Test all grid points
        x0s = X;

        for j = 1:length(fs)
                f = fs{j};
                fx = f(X(:,1), X(:,2));
            for i = 1:length(x0s)
                x0 = x0s(i,:);
                delta = diracDiscr(x0, x, mom_cond, 0, H);
                integral = delta'*H_global*fx;
                err = abs(integral - f(x0(1), x0(2)));
                testCase.verifyLessThan(err, 1e-12);
            end
        end
    end
end

function testAllRandom2D(testCase)

    orders = [2, 4, 6];
    mom_conds = orders;
    rng(0xDABBAD00) % Set seed. Jabba-dabba-doooo!

    for o = 1:length(orders)
        order = orders(o);
        mom_cond = mom_conds(o);
        [xlims, ylims, x, X, h, H, fs] = setup2D(order, mom_cond);
        H_global = kron(H{1}, H{2});

        xl = xlims{1};
        xr = xlims{2};
        yl = ylims{1};
        yr = ylims{2};

        % Test random points, even outside grid
        Npoints = 100;
        x0s = [(xl-3*h{1}) + (xr-xl+6*h{1})*rand(Npoints,1), ...
               (yl-3*h{2}) + (yr-yl+6*h{2})*rand(Npoints,1) ];

        for j = 1:length(fs)
                f = fs{j};
                fx = f(X(:,1), X(:,2));
            for i = 1:length(x0s)
                x0 = x0s(i,:);
                delta = diracDiscr(x0, x, mom_cond, 0, H);
                integral = delta'*H_global*fx;

                % Integral should be 0 if point is outside grid
                if x0(1) < xl || x0(1) > xr || x0(2) < yl || x0(2) > yr
                    err = abs(integral - 0);
                else
                    err = abs(integral - f(x0(1), x0(2)));
                end
                testCase.verifyLessThan(err, 1e-12);
            end
        end
    end
end

%=============== 3D tests ==============================
function testAllGP3D(testCase)

    orders = [2, 4, 6];
    mom_conds = orders;

    for o = 1:length(orders)
        order = orders(o);
        mom_cond = mom_conds(o);
        [~, ~, ~, x, X, ~, H, fs] = setup3D(order, mom_cond);
        H_global = kron(kron(H{1}, H{2}), H{3});

        % Test all grid points
        x0s = X;

        for j = 1:length(fs)
                f = fs{j};
                fx = f(X(:,1), X(:,2), X(:,3));
            for i = 1:length(x0s)
                x0 = x0s(i,:);
                delta = diracDiscr(x0, x, mom_cond, 0, H);
                integral = delta'*H_global*fx;
                err = abs(integral - f(x0(1), x0(2), x0(3)));
                testCase.verifyLessThan(err, 1e-12);
            end
        end
    end
end

function testAllRandom3D(testCase)

    orders = [2, 4, 6];
    mom_conds = orders;
    rng(0xDABBAD00) % Set seed. Jabba-dabba-doooo!

    for o = 1:length(orders)
        order = orders(o);
        mom_cond = mom_conds(o);
        [xlims, ylims, zlims, x, X, h, H, fs] = setup3D(order, mom_cond);
        H_global = kron(kron(H{1}, H{2}), H{3});

        xl = xlims{1};
        xr = xlims{2};
        yl = ylims{1};
        yr = ylims{2};
        zl = zlims{1};
        zr = zlims{2};

        % Test random points, even outside grid
        Npoints = 200;
        x0s = [(xl-3*h{1}) + (xr-xl+6*h{1})*rand(Npoints,1), ...
               (yl-3*h{2}) + (yr-yl+6*h{2})*rand(Npoints,1), ...
               (zl-3*h{3}) + (zr-zl+6*h{3})*rand(Npoints,1) ];

        for j = 1:length(fs)
                f = fs{j};
                fx = f(X(:,1), X(:,2), X(:,3));
            for i = 1:length(x0s)
                x0 = x0s(i,:);
                delta = diracDiscr(x0, x, mom_cond, 0, H);
                integral = delta'*H_global*fx;

                % Integral should be 0 if point is outside grid
                if x0(1) < xl || x0(1) > xr || x0(2) < yl || x0(2) > yr || x0(3) < zl || x0(3) > zr
                    err = abs(integral - 0);
                else
                    err = abs(integral - f(x0(1), x0(2), x0(3)));
                end
                testCase.verifyLessThan(err, 1e-12);
            end
        end
    end
end


% ======================================================
% ============== Setup functions =======================
% ======================================================
function [xl, xr, h, x, H, fs] = setup1D(order, mom_cond)

    % Grid
    xl = -3;
    xr = 900;
    L = xr-xl;
    m = 101;
    h = (xr-xl)/(m-1);
    g = grid.equidistant(m, {xl, xr});
    x = g.points();

    % Quadrature
    ops = sbp.D2Standard(m, {xl, xr}, order);
    H = ops.H;

    % Moment conditions
    fs = cell(mom_cond,1);
    for p = 0:mom_cond-1
        fs{p+1} = @(x) (x/L).^p;
    end

end

function [xlims, ylims, x, X, h, H, fs] = setup2D(order, mom_cond)

    % Grid
    xlims = {-3, 20};
    ylims = {-11,5};
    Lx = xlims{2} - xlims{1};
    Ly = ylims{2} - ylims{1};

    m = [15, 16];
    g = grid.equidistant(m, xlims, ylims);
    X = g.points();
    x = g.x;

    % Quadrature
    opsx = sbp.D2Standard(m(1), xlims, order);
    opsy = sbp.D2Standard(m(2), ylims, order);
    Hx = opsx.H;
    Hy = opsy.H;
    H = {Hx, Hy};

    % Moment conditions
    fs = cell(mom_cond,1);
    for p = 0:mom_cond-1
        fs{p+1} = @(x,y) (x/Lx + y/Ly).^p;
    end

    % Grid spacing in interior
    mm = round(m/2);
    hx = x{1}(mm(1)+1) - x{1}(mm(1));
    hy = x{2}(mm(2)+1) - x{2}(mm(2));
    h = {hx, hy};

end

function [xlims, ylims, zlims, x, X, h, H, fs] = setup3D(order, mom_cond)

    % Grid
    xlims = {-3, 20};
    ylims = {-11,5};
    zlims = {2,4};
    Lx = xlims{2} - xlims{1};
    Ly = ylims{2} - ylims{1};
    Lz = zlims{2} - zlims{1};

    m = [13, 14, 15];
    g = grid.equidistant(m, xlims, ylims, zlims);
    X = g.points();
    x = g.x;

    % Quadrature
    opsx = sbp.D2Standard(m(1), xlims, order);
    opsy = sbp.D2Standard(m(2), ylims, order);
    opsz = sbp.D2Standard(m(3), zlims, order);
    Hx = opsx.H;
    Hy = opsy.H;
    Hz = opsz.H;
    H = {Hx, Hy, Hz};

    % Moment conditions
    fs = cell(mom_cond,1);
    for p = 0:mom_cond-1
        fs{p+1} = @(x,y,z) (x/Lx + y/Ly + z/Lz).^p;
    end

    % Grid spacing in interior
    mm = round(m/2);
    hx = x{1}(mm(1)+1) - x{1}(mm(1));
    hy = x{2}(mm(2)+1) - x{2}(mm(2));
    hz = x{3}(mm(3)+1) - x{3}(mm(3));
    h = {hx, hy, hz};

end