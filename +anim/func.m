% Animates a function F(x,t) for x in lim from t to tend
function func(F, xrange, yrange, T, timeMod)
    default_arg('timeMod', 1);
    x = linspace(xrange(1),xrange(2), 200);

    fig_handle = figure;
    plot_handle = plot(x,F(x,T(1)));
    xlim(xrange)
    ylim(yrange)
    axis_handle = gca;


    function t = G(t)
        set(plot_handle,'YData',F(x,t))
        title(axis_handle,sprintf('T=%.3f',t));
        drawnow
    end

    anim.animate(@G, T(1), T(2), timeMod);
end





