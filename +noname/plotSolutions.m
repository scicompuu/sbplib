function plotSolutions(filename, figename_prefix, plot_type)
    default_arg('figename_prefix',[]);
    default_arg('plot_type','plot')

    save_figures = isempty(figename_prefix);

    sf = SolutionFile(filename);

    for i = 1:length(sf.keys)
        key = sf.keys{i};
        entry = sf.get(key);

        method  = key.method;
        order   = key.order;
        m       = key.m;
        T       = key.T;
        repr       = entry.repr;
        runtime = entry.runtime;
        k       = entry.k;

        discr = entry.discrHand(m,order);

        [update, hand] = discr.setupPlot(plot_type);
        update(repr);
    end

end