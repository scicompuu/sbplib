function [update_data,figure_handle] = setup_2d_plot(x,y,z_lim,zfun)
    default_arg('zfun',@(z)z);

    Z = zeros(length(y),length(x));
    figure_handle = figure;
    plot_handle = surf(x,y,Z);
    plot_handle.LineStyle = 'none';
    axis_handle = gca;

    xlabel('x')
    ylabel('y')
    xlim([x(1) x(end)]);
    ylim([y(1) y(end)]);
    zlim(z_lim);
    caxis(z_lim);
    % axis vis3d
    % colormap(parula(256))
    % colorbar

    function update(t,z)
        Z = zfun(z);
        % Z = reshape(zfun(z),length(x),length(y));
        if ishandle(plot_handle) && ishandle(axis_handle)
            set(plot_handle,'ZData',Z)
            title(axis_handle,sprintf('T=%.3f',t));
            drawnow
        end
    end
    update_data = @update;
end
