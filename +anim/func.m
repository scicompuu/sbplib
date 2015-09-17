% Animates a function F(x,t) for x in lim from t to tend
function func(xrange,yrange,t,tend, F)
    x = linspace(xrange(1),xrange(2), 200);

    fig_handle = figure;
    plot_handle = plot(x,F(x,t));
    xlim(xrange)
    ylim(yrange)
    axis_handle = gca;


    function t = G(t)
        set(plot_handle,'YData',F(x,t))
        title(axis_handle,sprintf('T=%.3f',t));
        drawnow
    end

    anim.animate(@G, t, tend)
end





