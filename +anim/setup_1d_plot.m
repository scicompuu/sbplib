% Creates a plot and provides a function to update the data in it.
%   x     - Vector of x-values to plot for.
%   y_lim - 1x2 vector containing the y limits of the plot.
%   yfun  - Function or a cell array of functions of y data vectors
%           that should be plotted. The output of each function
%           will be plotted to the same axis.
%
%   update_data(t,varargin) - Function to update plot data. All varargin will
%                             be passed to functions in yfun.
%   plot_handles            - Array of plot_handles. One for each yfun.
%   axis_handle             - Handle to the axis plotted to.
function [update_data, plot_handles, axis_handle] = setup_1d_plot(x,yfun,show_title)
    default_arg('yfun',{@(y)y});
    default_arg('show_title', true)

    if isa(yfun,'function_handle')
        yfun = {yfun};
    end

    figure_handle = gcf;
    plot_handles(1) = plot(x,0*x);
    for i = 2:length(yfun)
        plot_handles(i) = line(x,0*x);
    end

    axis_handle = gca;

    xlim([x(1) x(end)]);

    function update(t,varargin)
        if ishandle(figure_handle) && ishandle(axis_handle)
            for i = 1:length(yfun)
                set(plot_handles(i),'YData',yfun{i}(varargin{:}));
            end

            if show_title
                title(axis_handle,sprintf('T=%.3f',t));
            end
        end
    end
    update_data = @update;
end
