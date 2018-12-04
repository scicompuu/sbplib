% calculate a N-D mononomial with powers k in points x:
%  z = x(:,1).^k(1) * x(:,2).^k(2) * ...
function z = mononomial(x, k)
    assert(size(x,2) == length(k), 'k must have the same length as the width of x');

    if any(k < 0)
        z = x(:,1)*0;
        return
    end

    denom = prod(factorial(k));

    for i = 1:length(k)
        x(:,i) = x(:,i).^k(i);
    end
    z = prod(x,2)/denom;
end
