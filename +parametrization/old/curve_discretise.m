% Discretises the curve g with the smallest number of points such that all segments
% are shorter than h. If do_plot is true the points of the discretisation and
% the normals of the curve in those points are plotted.
%
%   [t,p,d] = curve_discretise(g,h,do_plot)
%
%   t is a vector of input values to g.
%   p is a cector of points.
%   d are the length of the segments.
function [t,p,d] = curve_discretise(g,h,do_plot)
    default_arg('do_plot',false)

    n = 10;

    [t,p,d] = curve_discretise_n(g,n);

    % ni = 0;
    while any(d>h)
        [t,p,d] = curve_discretise_n(g,n);
        n = ceil(n*d(1)/h);
        % ni = ni+1;
    end

    % nj = 0;
    while all(d<h)
        [t,p,d] = curve_discretise_n(g,n);
        n = n-1;
        % nj = nj+1;
    end
    [t,p,d] = curve_discretise_n(g,n+1);

    % fprintf('ni = %d, nj = %d\n',ni,nj);

    if do_plot
        fprintf('n:%d  max: %f min: %f\n', n, max(d),min(d));
        p = parametrization.map_curve(g,t);
        figure
        show(g,t,h);
    end

end

function [t,p,d] = curve_discretise_n(g,n)
    t = linspace(0,1,n);
    t = equalize_d(g,t);
    d = D(g,t);
    p = parametrization.map_curve(g,t);
end

function d = D(g,t)
    p = parametrization.map_curve(g,t);

    d = zeros(1,length(t)-1);
    for i = 1:length(d)
        d(i) = norm(p(:,i) - p(:,i+1));
    end
end

function t = equalize_d(g,t)
    d = D(g,t);
    v = d-mean(d);
    while any(abs(v)>0.01*mean(d))
        dt = t(2:end)-t(1:end-1);
        t(2:end) = t(2:end) - cumsum(dt.*v./d);

        t = t/t(end);
        d = D(g,t);
        v = d-mean(d);
    end
end


function show(g,t,hh)
    p = parametrization.map_curve(g,t);



    h = parametrization.plot_curve(g);
    h.LineWidth = 2;
    axis equal
    hold on
    h = plot(p(1,:),p(2,:),'.');
    h.Color = [0.8500 0.3250 0.0980];
    h.MarkerSize = 24;
    hold off

    n = parametrization.curve_normals(g,t);
    hold on
    for  i = 1:length(t)
        p0 = p(:,i);
        p1 = p0 + hh*n(:,i);
        l = [p0, p1];
        h = plot(l(1,:),l(2,:));
        h.Color = [0.8500 0.3250 0.0980];
    end

end