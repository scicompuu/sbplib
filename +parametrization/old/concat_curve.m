% Concatenate two curves g1 and g2 intop
%   g = concat_curve(g1,g2)
function g = concat_curve(g1,g2)
    function v = g_fun(t)
        if t < 1/2
            v = g1(2*t);
        else
            v = g2(2*t-1);
        end
    end
    g = @g_fun;
end