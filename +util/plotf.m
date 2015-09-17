function plotf(f,a,b,opt)
    if ~exist('opt')
        opt = '';
    end

    x = linspace(a,b);
    for i = 1:length(x)
        y(i) = f(x(i));
    end
    plot(x,y,opt)
end