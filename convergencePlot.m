function hand = convergencePlot(orders, h, e)
    N = length(orders);

    fh = figure();
    ah = axes();
    ah.XScale = 'log';
    ah.YScale = 'log';
    hold on
    ph = {};
    phc = {};
    legends = {};
    for i = 1:N
        ph{i} = loglog(h{i}, e{i});
        phc{i} = plotConvergenceFit(orders{i}, h{i}, e{i});

        ph{i}.LineStyle = 'none';
        ph{i}.Marker = Color.solidMarkers{i};
        ph{i}.MarkerSize = 12;
        ph{i}.Color = Color.colors{i};
        ph{i}.MarkerFaceColor = Color.colors{i};

        legends{i} = sprintf('$o = %d$', orders{i});
    end
    hold off

    lh = legend([ph{:}], legends);
    lh.Interpreter = 'latex';
    lh.Location = 'SouthEast';

    for i = 1:N
        uistack(phc{i}, 'bottom');
    end

    xlabel('$h$', 'interpreter', 'latex')
    ylabel('Error', 'interpreter', 'latex')

    % xlim([0.7e-2, 1e-1])
    % ylim([3e-5, 4])

    grid on

    ah = gca();
    ah.TickLabelInterpreter = 'latex';
    setFontSize(fh);

    % if savePngs
    %     savepng(fh, 'fig/conv/conv',600)
    % end

    hand = struct();
    hand.fig = fh;
    hand.data = ph;
    hand.fits = phc;
    hand.legend = lh;
end
