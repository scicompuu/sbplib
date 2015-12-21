function place_label(pt,str,horzAl,vertAl)
    default_arg('horzAl','center');
    default_arg('vertAl', 'middle');

    x = pt(1);
    y = pt(2);
    h = text(x,y,str);
    h.HorizontalAlignment = horzAl;
    h.VerticalAlignment = vertAl;
end