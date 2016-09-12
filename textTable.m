
% data           -- cell array of numbers
% leftColstrings -- cell array of strings, for left column
% topRowStrings  -- cell array of strings, for top row
% dataFormat     -- (optional) format specifier, e.g. '%.2f' 
function textTable(data, leftColStrings, topRowStrings, dataFormat)

    default_arg('dataFormat','%.2f')
    
    nRows = length(leftColStrings);
    nCols = length(topRowStrings);
    [m,n] = size(data);
    
    if(m ~= nRows || n ~=nCols)
        error('Data dimensions must match labels');
    end

    % Find column widths
    headerLength = 0;
    for i = 1:nCols
        headerLength = max(headerLength, length(topRowStrings{i} ));
    end
 
    dataLength = 0;
    for i = 1:nRows
        for j = 1:nCols
            temp = length(sprintf(dataFormat, data{i,j}));
            dataLength = max(dataLength, temp);
        end
    end
    dataLength = length(sprintf(dataFormat, data{1,1}));
    
    colWidth = max(headerLength,dataLength);

    % Print headers
    fprintf(' %*s |',colWidth,'')
    for i = 1:nCols
        fprintf(' %-*s |', colWidth, topRowStrings{i});
    end
    fprintf('\n');

    % Print divider
    m_dev = repmat('-',1,colWidth);
    column_dev = repmat('-',1,colWidth);
    fprintf('-%s-+',m_dev);
    for i = 1:nCols
        fprintf('-%s-+', column_dev);
    end
    fprintf('\n');
    

    % Print each row
    dataFormat = ['%*' dataFormat(2:end)];
    for i = 1:nRows
        fprintf(' %*s |',colWidth,leftColStrings{i});
        for j = 1:nCols
            fprintf([' ' dataFormat ' |'], colWidth, data{i,j});
        end
        fprintf('\n');
    end

    fprintf('\n');

end