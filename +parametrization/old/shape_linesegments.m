% Converts a shape into a cell array of linesegments shorter than h.
function l = shape_linesegments(s,h)
    l = {};

    for i = 1:length(s)
        t = parametrization.curve_discretise(s{i},h);
        l = [l, parametrization.curve_linesegments(s{i},t)];
    end
end