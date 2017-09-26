% Differentiates a symbolic function like diff does, but keeps the function as a symfun
function g = diffSymfun(f, varargin)
    assertType(f, 'symfun');

    args = argnames(f);
    g = symfun(diff(f,varargin{:}), args);
end
