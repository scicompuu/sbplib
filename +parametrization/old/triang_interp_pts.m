% Creates a transfinite interpolation from connecting the four points wiht straight lines.
function [S, g1, g2, g3] = triang_interp_pts(p1,p2,p3)
    if size(p1) ~= [2 1]
        error('p1 is strange!');
    end

    g1 = @(t)(p1 + t*(p2-p1));
    g2 = @(t)(p2 + t*(p3-p2));
    g3 = @(t)(p3 + t*(p1-p3));

    S = parametrization.triang_interp(g1,g2,g3);
end
