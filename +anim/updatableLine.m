function [update_data, line_handle] = updatableLine()
    figure_handle = gcf;
    axis_handle = gca;

    line_handle = line(0,0);

    function update(x,y)
        if ishandle(figure_handle) && ishandle(axis_handle)
            set(line_handle,'XData',x);
            set(line_handle,'YData',y);
            drawnow
        end
    end
    update_data = @update;
end