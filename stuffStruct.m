function s = stuffStruct(varargin)
    s = struct();

    for i = 1:nargin
        assert(~isempty(inputname(i)), 'All inputs must be variables.');
        s.(inputname(i)) = varargin{i};
    end
end
