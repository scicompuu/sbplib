% spyh mimics the built in spy but returns a handle which allows modifying the plot.
function h = spyh(A)
    [n,m] = size(A);
    [I,J] = find(A);

    if ~isempty(J)
        h = plot(J,I);
        h.LineStyle = 'none';
        h.Marker = '.';
        h.MarkerSize = 14;
    else
        h = [];
    end

    a = gca;
    xlim([0 m+1]);
    ylim([0 n+1]);
    axis square
    a.YDir = 'reverse';
end