function [update_data, plot_handles] = setup_time_quantity_plot(yfun)
    default_arg('yfun',@(y)y);

    if isa(yfun,'function_handle')
        yfun = {yfun};
    end

    t = [];
    for i = 1:length(yfun)
        plot_handles(i) = line(0,0);
        plot_handles(i).XData = [];
        plot_handles(i).YData = [];
        quantities{i} = [];
    end

    axis_handle = gca;
    legend()


    function update(t_now,varargin)
        if ishandle(axis_handle)
            t = [t t_now];
            for j = 1:length(yfun)
                quantities{j} = [quantities{j} yfun{j}(varargin{:})];
                plot_handles(j).XData = t;
                plot_handles(j).YData = quantities{j};
            end

            if t(end) > t(1)
                xlim([t(1) 1.1*t(end)]);
            end
        end
    end
    update_data = @update;
end
