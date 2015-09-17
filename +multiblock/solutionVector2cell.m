function v_cell = solutionVector2cell(schms,v,components)
    if nargin == 2
        components = 1;
    end

    N = length(schms);
    v_cell = cell(N,1);

    i_start = 1;
    for i =1:N
        i_end = i_start + schms{i}.size()*components - 1;
        v_cell{i} = v(i_start:i_end);
        i_start = i_end+1;
    end
end