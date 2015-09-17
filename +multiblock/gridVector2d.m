function [x,y] = gridVector2d(schms)
    n = length(schms);
    x = cell(n,1);
    y = cell(n,1);

    for i = 1:n
        x{i} = schms{i}.x;
        y{i} = schms{i}.y;
    end
end