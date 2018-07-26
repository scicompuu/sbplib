function [update_data, plot_handles] = setup_time_quantity_plot(yfun)
    default_arg('yfun',@(y)y);

    if isa(yfun,'function_handle')
        yfun = {yfun};
    end

    t = [];
    for i = 1:length(yfun)
        plot_handles(i) = animatedline();
    end

    axis_handle = gca;

    function update(t_now,varargin)
        if ishandle(axis_handle)
            % t = [t t_now];
            for j = 1:length(yfun)
                addpoints(plot_handles(j),t_now,full(yfun{j}(varargin{:})));
            end

            [t,~] = getpoints(plot_handles(1));
            if t(1) < t(end)
                xlim(axis_handle, [t(1) t(end)]);
            end
        end
    end
    update_data = @update;
end
