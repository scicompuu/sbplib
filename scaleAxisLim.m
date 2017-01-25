function scaleAxisLim(ah, s)
    x = ah.XLim;
    y = ah.YLim;

    xl = x(2) - x(1);
    yl = y(2) - y(1);

    dx = xl*(s-1);
    dy = yl*(s-1);

    ah.XLim = x + [-dx/2 dx/2];
    ah.YLim = y + [-dy/2 dy/2];
end