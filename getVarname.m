function names = getVarname(varargin)
    names = cell(size(varargin));

    for i = 1:numel(varargin)
        names{i} = inputname(i);
    end
end