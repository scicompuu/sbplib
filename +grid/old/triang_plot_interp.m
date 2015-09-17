% Plots a transfinite interpolation in x,y space using nu and nv curves along u and v axes.






% Plots a interp of a triangle where one the interpolation is from a square
% with one side collapsed to
function h = triang_plot_interp_kindaworking(S,n)
    u = linspace(0,1,n);
    v = linspace(0,1,n);

    m = 100;
    m = 20;

    Xl_curves = cell(n,1);
    Xr_curves = cell(n,1);
    Y_curves = cell(n,1);


    function u = wierdness(v,d,N)
        if N == 0
            u = 0;
        else
            u = N*d./(1-v);
        end
    end


    %Y curves
    t = linspace(0,1,m);
    for i = 1:n
        x = []; y = [];
        for j = 1:length(t)
            [x(j),y(j)] = S(t(j),v(i));
        end
        Y_curves{i} = [x', y'];
    end


    % Right and left X curves
    t = linspace(0,1,m);
    d = u(2);
    for i = 1:n
        xl = []; yl = [];
        xr = []; yr = [];
        N = i-1;
        t = linspace(0,1-N*d,m);
        for j = 1:length(t)
            w = wierdness(t(j),d,N);
            [xr(j),yr(j)] = S(w,t(j));
            [xl(j),yl(j)] = S(1-w,t(j));
        end
        Xl_curves{i} = [xl', yl'];
        Xr_curves{i} = [xr', yr'];
    end

    for i = 1:n-1
        line(Xl_curves{i}(:,1),Xl_curves{i}(:,2))
        line(Xr_curves{i}(:,1),Xr_curves{i}(:,2))
        line(Y_curves{i}(:,1),Y_curves{i}(:,2))
    end
end




function h = triang_plot_interp_nonworking(S,n)

    u = linspace(0,1,n);
    v = linspace(0,1,n);

    m = 100;

    X_curves = cell(n-1,1);
    Y_curves = cell(n-1,1);
    K_curves = cell(n-1,1);


    t = linspace(0,1,m);
    for i = 1:n-1
        x = []; y = [];
        for j = find(t+u(i) <= 1)
            [x(j),y(j)] = S(u(i),t(j));
        end
        X_curves{i} = [x', y'];
    end

    for i = 1:n-1
        x = []; y = [];
        for j = find(t+v(i) <= 1)
            [x(j),y(j)] = S(t(j),v(i));
        end
        Y_curves{i} = [x', y'];
    end

    for i = 2:n
        x = []; y = [];
        for j = find(t<u(i))
            [x(j),y(j)] = S(t(j), u(i)-t(j));
        end
        K_curves{i-1} = [x', y'];
    end

    for i = 1:n-1
        line(X_curves{i}(:,1),X_curves{i}(:,2))
        line(Y_curves{i}(:,1),Y_curves{i}(:,2))
        line(K_curves{i}(:,1),K_curves{i}(:,2))
    end

    h = -1;
    % h = plot(X_curves{:},Y_curves{:},K_curves{:});
end
