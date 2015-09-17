function [p, button] = prompt_point(s,varargin)
    default_arg('s',[])

    set(gcf,'Pointer','crosshair')

    if ~isempty(s)
        fprintf(s,varargin{:});
    end

    a = gca;

    function get_point(src,event)
        cp = a.CurrentPoint;
        p = cp(1,1:2)';
        a.ButtonDownFcn = [];
    end

    a.ButtonDownFcn = @get_point;
    waitfor(a,'ButtonDownFcn', [])

    set(gcf,'Pointer','arrow')

end