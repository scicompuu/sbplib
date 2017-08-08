function S = spzeros(varargin)
    switch length(varargin)
        case 2
            S = sparse(varargin{1}, varargin{2});
        case 1
            v = varargin{1};
            switch length(v)
                case 1
                    S = sparse(v,v);
                case 2
                    S = sparse(v(1), v(2));
                otherwise
                    error('Input must be either one integer, two integers or a vector with two integers');
            end
        otherwise
            error('Too many input arguments.');
    end
end
