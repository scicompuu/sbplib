function x = gridVector1d(schms)
    n = length(schms);
    x = cell(n,1);

    for i = 1:n
        x{i} = schms{i}.x;
    end
end