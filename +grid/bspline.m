% Calculates a D dimensional p-order bspline at t with knots T and control points P.
%  T = [t0 t1 t2 ... tm] is a 1 x (m+1) vector with non-decresing elements and t0 = 0 tm = 1.
%  P = [P0 P1 P2 ... Pn] is a D x (n+1) matrix.

% knots p+1 to m-p-1 are the internal knots

% Implemented from: http://mathworld.wolfram.com/B-Spline.html
function C = bspline(t,p,P,T)
    m = length(T) - 1;
    n = size(P,2) - 1;
    D = size(P,1);

    assert(p == m - n - 1);

    C = zeros(D,length(t));

    for i = 0:n
        for k = 1:D
            C(k,:) = C(k,:) + P(k,1+i)*B(i,p,t,T);
        end
    end

    % Curve not defined for t = 1 ? Ugly fix:
    I = find(t == 1);
    C(:,I) = repmat(P(:,end),[1,length(I)]);
end

function o = B(i, j, t, T)
    if j == 0
        o = T(1+i) <= t & t < T(1+i+1);
        return
    end

    if T(1+i+j)-T(1+i) ~= 0
        a = (t-T(1+i))/(T(1+i+j)-T(1+i));
    else
        a = t*0;
    end

    if T(1+i+j+1)-T(1+i+1) ~= 0
        b = (T(1+i+j+1)-t)/(T(1+i+j+1)-T(1+i+1));
    else
        b = t*0;
    end

    o = a.*B(i, j-1, t, T) + b.*B(i+1, j-1, t, T);
end