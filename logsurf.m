function [srfHandle, cbHandle] = logsurf(X,Y,Z, lim)
    absLogZ = log10(abs(Z));
    srfHandle = surf(X,Y,absLogZ);

    cbHandle = colorbar();
    colormap(hot(256));
    ah = gca();
    ah.CLim = lim;

    oldTickLabels = cbHandle.TickLabels;

    newTickLabels = {};

    for i = 1:length(oldTickLabels)
        newTickLabels{i} = sprintf('10^{%s}',oldTickLabels{i});
    end

    cbHandle.TickLabels = newTickLabels;
end