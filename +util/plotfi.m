% runs f for each value in v and plots the result.
function plotfv(f,v,opt)

    if ~exist('opt')
        opt = '';
    end

    for i = 1:length(v)
        y(i) = f(v(i));
    end
    plot(v,y,opt)
end