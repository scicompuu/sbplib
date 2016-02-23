function ineq = semiDefIneq(A)
    [m, sub] = minors(A);

    ineqsys = true;
    for i = 1:length(m)
        ineqsys = ineqsys & m(i) >= 0;
    end

    ineq = simplify(ineqsys);

    str = toString(ineq);
    fprintf('%s\n',strjoin(strsplit(str,' & '), '\n'));
end