function assertNumberOfArguments(fun, N)
    if nargin(fun) ~= N
        error('sbplib:assertNumberOfArguments:wrongNumberOfArguments', '"%s" must have %d, found %d', inputname(1), N, nargin(fun));
    end
end
