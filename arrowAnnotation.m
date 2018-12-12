% Draw an arrow from p1 to p2, with text attached
function [h] = arrowAnnotation(p1,p2)
    ah = gca;
    xl = ah.XLim(1);
    xr = ah.XLim(2);

    yl = ah.YLim(1);
    yr = ah.YLim(2);

    dx = xr - xl;
    dy = yr - yl;

    s = [
        ah.Position(1) + (p1(1) - xl)/dx*ah.Position(3),
        ah.Position(1) + (p2(1) - xl)/dx*ah.Position(3),
    ];
    t = [
        ah.Position(2) + (p1(2) - yl)/dy*ah.Position(4),
        ah.Position(2) + (p2(2) - yl)/dy*ah.Position(4),
    ];

    h = annotation('arrow', s, t);
end
