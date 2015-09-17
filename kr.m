function A = kr(varargin)
    n = nargin;

    A = 1;

    for i = 1:n
        A = kron(A,varargin{i});
    end
end