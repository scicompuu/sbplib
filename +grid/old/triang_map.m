% Creates a grid [X,Y] from the mapping function S at points in vectors u,v
function [X, Y] = traing_map(S,u,v)
    error('not done')
    if nargin == 2
        v = u;
    end

    if isscalar(u)
        u = linspace(0,1,u);
    end

    if isscalar(v)
        v = linspace(0,1,v);
    end

    nu = length(u);
    nv = length(v);

    X = zeros(nu,nv);
    Y = zeros(nu,nv);

    for i = 1:nu
        for j = 1:nv
            [x,y] = S(u(i),v(j));
            X(i,j) = x;
            Y(i,j) = y;
        end
    end
end
