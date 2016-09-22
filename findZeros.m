% findZeros looks for solutions to the equation f(x)==0 within
% the limits lim with a granularity of h.
% Returns a sorted list of unique solutions.
function z = findZeros(f, lim, h)
    n = ceil((lim(2)-lim(1))/h);
    z0 = linspace(lim(1), lim(2), n);

    z = zeros(1,n);

    for i = 1:n
        zt(i) = fzero(f, z0(i));
    end

    zt = sort(zt);

    z = [];
    for i = 1:n
        if zt(i)  < lim(1) || zt(i) > lim(2)
            continue
        end

        if ~isempty(z) && abs(z(end) - zt(i)) < 1e-6
            continue
        end

        z = [z zt(i)];
    end

    % z = unique(z);
end