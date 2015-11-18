% [x,h] = get_grid(a,b,m)
function [x,h] = get_grid(a,b,m)
    % TODO: allow the interval to be a vector
    x = linspace(a,b,m)';
    h = (b-a)/(m-1);
end