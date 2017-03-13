% Draws a line in a loglog plot with slope `slope` fitted to the data in `x`
% and `y`. `xlength` scales how much of the interval [x(1) x(end)] is coverd
% by the line. `offset` is a multiplicative offset to where the line is drawn
% relative to the data.
function hand = plotConvergenceFit(slope, x, y, xlength, offset)
    default_arg('xlength', 0.8)
    default_arg('offset', 1);

    % Optimise for log(y) = p*log(x) + q

    p = slope;

    logx = log(x);
    logy = log(y);

    N = length(logx);

    q = 1/N*sum(logy-p*logx);

    logxlength = xlength * abs(logx(end)-logx(1));
    logxends = (logx(1)+logx(end))/2 + [-logxlength/2, logxlength/2];

    xends = exp(logxends);
    yends = exp(q)*xends.^p;

    hand = line(xends, yends);
    hand.Color = Color.black;
    hand.LineStyle = '--';
    hand.LineWidth = 2;
end



% function hand = plotConvergenceFit(slope, pos, width)
%     x0 = pos(1);
%     y0 = pos(2);

%     x = [x0*10^-(width/2) x0*10^(width/2)];
%     y = x.^slope * x0^-slope * y0;

%     hand = line(x,y);
%     hand.Color = Color.black;
%     hand.LineStyle = '--';
%     hand.LineWidth = 2;
% end