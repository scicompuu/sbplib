% Saves a figure to an .eos file with corrected bounding box.
function saveeps(handle,filename)
    if length(filename) < 4 || ~strcmp(filename(end-3:end),'.eps')
        filename = [filename '.eps'];
    end

    handle_units = handle.Units;  % Save the current units to be able to restore

    % Copy size of figure in centimeters to a place where saveas will honor it
    handle.Units = 'centimeters';
    handle.PaperUnits = 'centimeters';
    handle.PaperPosition(3:4) = handle.Position(3:4);

    % Save as a bugged eps file.
    saveas(handle,filename,'epsc');

    handle.Units = handle_units; % Restore the old units

    % Correct the buggy eps file
    correct_stupid_matlab_bug(filename);
end

% Corrects the format of an eps file so that the bounding box is defined at the top of the file instead of
% at the bottom.
function correct_stupid_matlab_bug(filename)
    contents = fileread(filename);
    lines = strsplit(contents,'\n');

    % Find the line
    pagel = findPrefix(lines,'%%Pages:');
    boundl = findPrefix(lines,'%%BoundingBox:');


    if ~(length(pagel) == 2 && length(boundl) == 2)
        warning('Undexpected number of found lines: %d , %d\nNot correcting the file',pagel, boundl);
        return
    end

    if ~(strcmp(lines{pagel(1)},'%%Pages: (atend)') && strcmp(lines{boundl(1)},'%%BoundingBox: (atend)'))
        warning('Does the file really contain the error?\nNot correcting the file');
        return
    end

    % Overwrite the nasty lines with the nice ones.
    lines{pagel(1)} = lines{pagel(2)};
    lines{boundl(1)} = lines{boundl(2)};

    % Delete the duplicates
    lines(pagel(2)) = [];
    lines(boundl(2)) = [];


    %Rewrite the file
    contents = strjoin(lines,'\n');

    fh = fopen(filename,'w');
    fprintf(fh, '%s',contents);
    fclose(fh);
end

function I = findPrefix(lines, prefix)
    I = find(strncmp(lines,prefix,length(prefix)));
end
