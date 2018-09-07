function [p, button] = prompt_point(s, varargin)
    default_arg('s',[])

    set(gcf,'Pointer','crosshair')

    if ~isempty(s)
        fprintf(s, varargin{:});
    end

    fh = gcf();
    ah = gca();

    function get_point(src, event)
        cp = ah.CurrentPoint;
        p = cp(1,1:2)';
        fh.WindowButtonUpFcn = [];
    end

    fh.WindowButtonUpFcn = @get_point;
    waitfor(fh,'WindowButtonUpFcn', [])

    set(gcf,'Pointer','arrow')

end