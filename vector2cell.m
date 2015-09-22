% Splits column vector v into segments of length n and returns the result as a column cell array.
%   v  -- column vector to be split
%   n  -- number of elements in each part
%
%   cv -- cell array of vectors with lenght n(i)
function cv = vector2cell(v,n)
    cv = cell(length(n),1);

    ind = [0; cumsum(n)];
    for i = 1:length(n)
        ind_i = (ind(i)+1):ind(i+1);
        cv{i} = v(ind_i);
    end
end