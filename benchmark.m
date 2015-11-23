function benchmark(func, N)
    default_arg('N',100);

    tic
    profile on

    for i = 1:N
        func();
    end

    profile viewer
end