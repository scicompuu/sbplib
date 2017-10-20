% Takes a function and evaluates it on a grid to return a grid function in the
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

    if g.D ~= nargin(func)
        error('grid:evalOn:WrongNumberOfInputs', 'The number of inputs of the function must match the dimension of the domain.')
    end

    x = g.points();
    k = numberOfComponents(func, x);

    % Evaluate gf = func(x(:,1),x(:,2),...,x(:,dim));
    if(g.D == 1)
        gf = func(x);
    else
        eval_str = 'gf = func(x(:,1)';
        for i = 2:g.D
            eval_str = [eval_str, sprintf(',x(:,%d)',i)];
        end
        eval_str = [eval_str, ');'];
        eval(eval_str);
    end

    % Reorganize gf
    gf_temp = gf;
    gf = zeros(g.N*k, 1);
    for i = 1:k
        gf(i:k:end) = gf_temp((i-1)*g.N + 1 : i*g.N);
    end
end

% Find the number of vector components of func
function k = numberOfComponents(func, x)
    if size(x,1) ~= 0
        x0 = x(1,:);
    else
        x0 = num2cell(ones(1,size(x,2)));
    end

    dim = length(x0);
    % Evaluate f0 = func(x0(1),x0(2),...,x0(dim));
    if(dim == 1)
        f0 = func(x0);
    else
        eval_str = 'f0 = func(x0(1)';
        for i = 2:dim
            eval_str = [eval_str, sprintf(',x0(%d)',i)];
        end
        eval_str = [eval_str, ');'];
        eval(eval_str);
    end

    % k = number of components
    k = length(f0);

    if size(f0,2) ~= 1
        error('grid:evalOn:VectorValuedWrongDim', 'A vector valued function must be given as a column vector')
    end
end