% Create a cubic spline from the points (x,y) using periodic conditions.
%   g = curve_interp(x,y)
function g = curve_interp(x,y)
    default_arg('x',[0 2 2 1 1 0])
    default_arg('y',[0 0 2 2 1 1])
    % solve for xp and yp

    % x(t) = at^4 + bt^2+ct+d

    % a = xp1 -2x1 + 2x0 +  xp0
    % b = 3x1 -xp1 - 3x0 + 2xp0
    % c = xp0
    % d = x0

    assert(length(x) == length(y))
    n = length(x);
    A = spdiags(ones(n,1)*[2, 8, 2],-1:1,n,n);
    A(n,1) = 2;
    A(1,n) = 2;

    bx = zeros(n,1);
    for i = 2:n-1
        bx(i) = -6*x(i-1)+6*x(i+1);
    end
    bx(1) = -6*x(n)+6*x(2);
    bx(n) = -6*x(n-1)+6*x(1);

    by = zeros(n,1);
    for i = 2:n-1
        by(i) = -6*y(i-1)+6*y(i+1);
    end
    by(1) = -6*y(n)+6*y(2);
    by(n) = -6*y(n-1)+6*y(1);


    xp = A\bx;
    yp = A\by;

    x(end+1) = x(1);
    y(end+1) = y(1);

    xp(end+1) = xp(1);
    yp(end+1) = yp(1);

    function v = g_fun(t)
        t = mod(t,1);
        i = mod(floor(t*n),n) + 1;
        t = t * n -(i-1);
        X = (2*x(i)-2*x(i+1)+xp(i)+xp(i+1))*t.^3 + (-3*x(i)+3*x(i+1)-2*xp(i)-xp(i+1))*t.^2 + (xp(i))*t + x(i);
        Y = (2*y(i)-2*y(i+1)+yp(i)+yp(i+1))*t.^3 + (-3*y(i)+3*y(i+1)-2*yp(i)-yp(i+1))*t.^2 + (yp(i))*t + y(i);
        v = [X;Y];
    end

    g = @g_fun;
end



