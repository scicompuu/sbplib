% Discretises a shape such that points on the curves are no further than h appart.
function p = shape_discretise(s,h)
    p = [];
    for i = 1:length(s)
        [~,pt] = grid.curve_discretise(s{i},h);
        p = [p, pt];
    end
end