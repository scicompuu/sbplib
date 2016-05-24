function t = time(f)
    s = tic();

    f();

    if nargout == 1
        t = toc(s);
    else
        toc(s);
end