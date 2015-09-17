function [update_data,figure_handle] = setup_2d_plot(x,y,clim,zfun)
    default_arg('zfun',@(z)z);

    Z = zeros(length(x),length(y));
    figure_handle = figure;
    plot_handle = imagesc(x,y,Z);
    axis_handle = gca;

    xlabel('x')
    ylabel('y')
    xlim([x(1) x(end)]);
    ylim([y(1) y(end)]);
    caxis(clim);
    % axis vis3d
    colormap(parula(256))
    colorbar

    function update(t,z)
        Z = zfun(z);
        % Z = reshape(zfun(z),length(x),length(y));
        if ishandle(plot_handle) && ishandle(axis_handle)
            set(plot_handle,'CData',Z)
            title(axis_handle,sprintf('T=%.3f',t));
            drawnow
        end
    end
    update_data = @update;
end


% TODO
% This function is for squre grid.
% This function is for 2d image.
% Make one for 3d surface
% Make one for curvilinear grids using pcolor


function [update_data,figure_handle,plot_handles] = setup_1d_plot(x,y_lim,yfun)

    figure_handle = figure;
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



