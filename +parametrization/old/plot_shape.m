% Plot a shape using n points. Returns cell array of plot handles.
%   hs = plot_shape(s,n)
function hs = plot_shape(s,n)
    default_arg('n',100);

    hs = {};
    hold on
    for i = 1:length(s)
        hs{end+1} = parametrization.plot_curve(s{i},n);
    end
    hold off
end