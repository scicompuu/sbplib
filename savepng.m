% TODO
% Let print size in inches as input parameter
% Smaller boundingbox
function savepng(h, filename, dpi)
    default_arg('dpi', 300)

    handle_units = h.Units;  % Save the current units to be able to restore

    % Copy size of figure in centimeters to a place where saveas will honor it
    h.Units = 'centimeters';
    h.PaperUnits = 'centimeters';
    h.PaperPosition(3:4) = h.Position(3:4);

    % Save as a bugged eps file.
    print(h,filename,'-dpng',sprintf('-r%d',dpi));

    h.Units = handle_units; % Restore the old units



end
