% Draws a line in a loglog plot with slope `slope` and position `pos`
% and `width` orders of magnitude wide.
%    ex: pos = [1e-1 1e-4]

function hand = plotConvergenceFit(slope, pos, width)
    x0 = pos(1);
    y0 = pos(2);

    x = [x0*10^-(width/2) x0*10^(width/2)];
    y = x.^slope * x0^-slope * y0;

    hand = line(x,y);
    hand.Color = Color.black;
    hand.LineStyle = '--';
    hand.LineWidth = 2;
end