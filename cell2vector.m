% cell2vector accepts a column cell array of column vectors and returns a columnvector
% with the input concatenated. It also returns the number of elements in each vector.
%   cv -- column cell array with column vectors
%   v  -- vector of the concatenated vectors
%   n  -- number of elements in each vector before concatenation. Can be used with vector2cell().
function [v, n] = cell2vector(cv)
    v = [];
    n = zeros(length(cv),1);

    for i = 1:length(cv)
        n(i) = length(cv{i});
        v = [v; cv{i}];
    end
end


% IS THIS ONE REALLY NEEDED? JUST USE cell2sparse?

