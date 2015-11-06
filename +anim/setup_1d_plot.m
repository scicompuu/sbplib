function [update_data,plot_handles] = setup_1d_plot(x,y_lim,yfun)
    default_arg('yfun',{@(y)y});

    if isa(yfun,'function_handle')
        yfun = {yfun};
    end

    figure_handle = gcf;
    plot_handles(1) = plot(x,0*x);
    hold on
    for i = 2:length(yfun)
        plot_handles(i) = plot(x,0*x);
    end
    hold off

    axis_handle = gca;

    xlabel('x')
    ylabel('y')
    xlim([x(1) x(end)]);
    ylim(y_lim);

    function update(t,varargin)
        if ishandle(figure_handle) && ishandle(axis_handle)
            for i = 1:length(yfun)
                set(plot_handles(i),'YData',yfun{i}(varargin{:}));
            end
            title(axis_handle,sprintf('T=%.3f',t));
            drawnow
        end
    end
    update_data = @update;
end
