function savepng(h, filename, dpi)
    default_arg('dpi', 300)
    print(h,filename,'-dpng',sprintf('-r%d',dpi));
    % Let print size in inches as input parameter
    % Smaller boundingbox
end
