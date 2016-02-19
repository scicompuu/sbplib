% Takes a funciton and evaluates it on a grid to return a grid function in the
% form of a (n*k)x1 vector, where n is the number of grid points and k is the
% number of components of the function.
%      g -- Grid to evaluate on.
%   func -- Function to evaluate. May be a function handle or a constant. If
%           it is a vector value it has to be provided as a column vector,
function gf = evalOn(g, func)
    if ~isa(func, 'function_handle')
        % We should have a constant.
        if size(func,2) ~= 1
            error('grid:evalOn:VectorValuedWrongDim', 'A vector valued function must be given as a column vector')
        end

        gf = repmat(func,[g.N, 1]);
        return
    end
    % func should now be a function_handle

    % Get coordinates and convert to cell array for easier use as a parameter
    x = g.points();
    X = {};
    for i = 1:size(x, 2)
        X{i} = x(:,i);
    end

    % Find the number of components
    x0 = num2cell(x(1,:));
    f0 = func(x0{:});
    k = length(f0);

    if size(f0,2) ~= 1
        error('grid:evalOn:VectorValuedWrongDim', 'A vector valued function must be given as a column vector')
    end

    gf = func(X{:});
    gf = reshape(reshape(gf, [g.N, k])', [g.N*k, 1]); % Reorder so that componets sits together.
end