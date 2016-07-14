function ineq = semiDefIneq(A, verbose)
    default_arg('verbose', true);
    [m, sub] = minors(A, verbose);

    ineqsys = true;
    for i = 1:length(m)
        ineqsys = ineqsys & m(i) >= 0;
    end

    ineq = simplify(ineqsys);

    str = toString(ineq);
    fprintf('%s\n',strjoin(strsplit(str,' & '), '\n'));
end