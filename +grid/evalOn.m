% Takes a function and evaluates it on a grid to return a grid function in the
% form of a (n*k)x1 vector, where n is the number of grid points and k is the
% number of components of the function.
%      g -- Grid to evaluate on.
%   func -- Function to evaluate. May be a function handle or a constant. If
%           it is a vector value it has to be provided as a column vector,
function gf = evalOn(g, func)
    if ~isa(func, 'function_handle')
        % We should have a constant.
        assert(size(func,2) == 1,'grid:evalOn:VectorValuedWrongDim', 'A vector valued function must be given as a column vector');

        gf = repmat(func,[g.N, 1]);
        return
    end
    % func should now be a function_handle
    assert(g.D == nargin(func) || nargin(func) < 0,'grid:evalOn:WrongNumberOfInputs', 'The number of inputs of the function must match the dimension of the domain.')

    x = num2cell(g.points(),1);
    k = numberOfComponents(func);

    gf = func(x{:});
    gf = reorderComponents(gf, k);
end

% Find the number of vector components of func
function k = numberOfComponents(func)
    x0 = num2cell(ones(1,nargin(func)));
    f0 = func(x0{:});
    assert(size(f0,2) == 1, 'grid:evalOn:VectorValuedWrongDim', 'A vector valued function must be given as a column vector');
    k = length(f0);
end

% Reorder the components of the function to sit together
function gf = reorderComponents(a, k)
    N = length(a)/k;
    gf = zeros(N*k, 1);
    for i = 1:k
        gf(i:k:end) = a((i-1)*N + 1 : i*N);
    end
end
