% data           -- cell array of numbers
% leftColstrings -- cell array of strings, for left column
% topRowStrings  -- cell array of strings, for top row
% dataFormat     -- (optional) format specifier, e.g. '%.2f' 
function latexTable(data, leftColStrings, topRowStrings, dataFormat)

    default_arg('dataFormat','%8.2f')
    
    nRows = length(leftColStrings);
    nCols = length(topRowStrings);
    [m,n] = size(data);
    
    if(m ~= nRows || n ~=nCols)
        error('Data dimensions must match labels');
    end

    header = {
        '\begin{table}[H]'
        '\centering'
        ['\begin{tabular}{c' repmat('|c',1,nCols) '} &']
        headers(topRowStrings)
        '\hline'
    };

    footer = {
        '\end{tabular}'
        '\caption{DESCRIPTION.}'
        '\label{table:LABEL}'
        '\end{table}'
    };
    
    nlc = sprintf('\n');
    dataStr = '';
    for i = 1:nRows
        dataStr = [dataStr leftColStrings{i}]; %#ok<AGROW>
        for j = 1:nCols
            dataStr = [dataStr ' & ' sprintf(dataFormat,data{i,j}) ]; %#ok<AGROW>
        end
        if(i<nRows)
            dataStr = [dataStr '  \\ ' nlc]; %#ok<AGROW>
        end
    end

    header = strjoin(header', nlc);
    footer = strjoin(footer', nlc);

    table = strjoin({header, dataStr, footer}, nlc);
    fprintf('%s\n', table);
end


function s = headers(strings)
    s= [strings{1} ' '];
    nCols = length(strings);
    for i = 2:nCols
        s = [s '& ' strings{i} ' ']; %#ok<AGROW>
    end
    s = [s ' \\'];
end