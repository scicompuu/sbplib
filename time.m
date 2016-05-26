function t_out = time(f, n)
    default_arg('n',1);

    if n == 1
        t = timeOnce(f);
    else
        t = timeAvg(f, n);
    end

    if nargout == 1
        t_out = t;
    else
        fprintf('Elapsed time is %.6f seconds.\n', t)
    end
end

function t = timeOnce(f)
    s = tic();

    f();

    t = toc(s);
end


function t = timeAvg(f, n)
    s = tic();

    for i = 1:n
        f();
    end

    t = toc(s)/n;
end
