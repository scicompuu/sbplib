% Takes a variable number of points as inputs and plots them with a label to the current figure.
%   label_pt(p)
function label_pt(varargin)
    for i = 1:length(varargin)
        try
            placelabel(varargin{i},inputname(i));
        catch e
            error('Could not place label for input %d, ''%s''. Not a valid point',i,inputname(i))
        end
    end
end

function placelabel(pt,str)
    x = pt(1);
    y = pt(2);
    h = line(x,y);
    h.Marker = '.';
    h = text(x,y,str);
    h.HorizontalAlignment = 'center';
    h.VerticalAlignment = 'bottom';
end