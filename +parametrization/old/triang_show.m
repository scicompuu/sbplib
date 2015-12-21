% Show a grid as a matlab figure.
function triang_show(S,n)

    ih = ishold();

    hold on
    h_grid = grid.triang_plot_interp(S,n);
    h_bord = grid.triang_plot_interp(S,2);

    set(h_grid,'Color',[0 0.4470 0.7410]);
    set(h_bord,'Color',[0.8500 0.3250 0.0980]);
    set(h_bord,'LineWidth',2);

    % axis auto
    % axis equal
    % axis square

    if ih
        hold on
    else
        hold off
    end
end
