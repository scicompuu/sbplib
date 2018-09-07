function [C,h] = prompt_bezier(s,varargin)
    default_arg('s',[])
    if ~isempty(s)
        fprintf(s,varargin{:});
    end

    fprintf('# Draw bezier curve\n')
    a = draw.prompt_point('Enter starting point\n');
    p = draw.point(a);
    p.Color = Color.green;
    p.MarkerSize = 24;
    b = draw.prompt_point('Enter stopping point\n');
    p = draw.point(b);
    p.Color = Color.red;
    p.MarkerSize = 24;
    c1 = draw.prompt_point('Enter control point 1\n');
    p = draw.point(c1);
    p.Color = Color.yellow;
    p.MarkerSize = 16;
    c2 = draw.prompt_point('Enter control point 2\n');
    p = draw.point(c2);
    p.Color = Color.yellow;
    p.MarkerSize = 16;

    C = parametrization.Curve.bezier(a,c1,c2,b);
    fprintf('C = parametrization.Curve.bezier([%.3g; %.3g],[%.3g; %.3g],[%.3g; %.3g],[%.3g; %.3g]);\n',a,c1,c2,b)
    h = C.plot();
    uistack(h,'bottom');

end