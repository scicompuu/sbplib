% Get the index of the points p within the tall array of points ps
function [I, ok] = pointIndex(p, ps)
    [ok, I] = ismember(p,  ps, 'rows');
end
