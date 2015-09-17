function [C,h] = prompt_line(s,varargin)
    default_arg('s',[])
    if ~isempty(s)
        fprintf(s,varargin{:});
    end

    a = draw.prompt_point('Enter starting point\n');
    p = draw.point(a);
    p.Color = Color.green;
    p.MarkerSize = 24;
    b = draw.prompt_point('Enter stopping point\n');
    p = draw.point(b);
    p.Color = Color.red;
    p.MarkerSize = 24;

    C = grid.Curve.line(a,b);
    h = C.plot();
    uistack(h,'bottom');
end