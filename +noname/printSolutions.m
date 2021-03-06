function printSolutions(filename)
    sf = SolutionFile(filename);

    name  = {};
    order   = [];
    m       = [];
    T       = [];
    t       = [];
    runtime = [];
    k       = [];

    for i = 1:length(sf.keys)
        key = sf.keys{i};
        entry = sf.get(key);

        name    = [name    key.name];
        order   = [order   key.order];
        m       = [m       key.m];
        T       = [T       key.T];
        t       = [t       entry.repr.t];
        runtime = [runtime entry.runtime];
        k       = [k       entry.k];
    end

    nameW    = findFieldWidth('%s', name);
    orderW   = findFieldWidth('%d', order);
    mW       = findFieldWidth('%d', m);
    TW       = findFieldWidth('%d', T);
    tW       = findFieldWidth('%.3e', t);
    runtimeW = findFieldWidth('%.3f', runtime);
    kW       = findFieldWidth('%.4f', k);

    for i = 1:length(sf.keys)
        fprintf('[%*s: o=%-*d, m=%-*d, T=%-*d]: t=%-*.3e, runtime=%*.3f, k=%*.4f\n', nameW, name{i}, orderW,order(i),mW,m(i),TW,T(i), tW, t(i), runtimeW,runtime(i), kW, k(i));
    end

end