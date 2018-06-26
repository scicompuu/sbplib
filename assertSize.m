% Assert that array A has the size s.
function assertSize(A,varargin)
    if length(varargin) == 1
        s = varargin{1};
        errmsg = sprintf('Expected %s to have size %s, got: %s',inputname(1), toString(s), toString(size(A)));
        assert(all(size(A) == s), errmsg);
    elseif length(varargin) == 2
        dim = varargin{1};
        s = varargin{2};

        errmsg = sprintf('Expected %s to have size %d along dimension %d, got: %d',inputname(1), s, dim, size(A,dim));
        assert(size(A,dim) == s, errmsg);
    else
        error('Expected 2 or 3 arguments to assertSize()');
    end
end
