function l = point(p)
    l = line(p(1,:),p(2,:));
    l.LineStyle = 'none';
    l.Marker = '.';
end
