
% data -- cell array of strings for each cell
function textTable(data,divCol, divColWeight, divRow, divRowWeight)

    % Find column widths
    colWidth = {};
    for i = 1:size(data, 2)
        for j = 1:size(data, 1)
            c(j) = length(data{j,i});
        end
        colWidth{i} = max(c);
    end

    error('not done')

    eW = 0;
    qW = 0;
    for i = 1:length(orders)
        log_e{i} = log10(e{i});
        eW = max(eW, findFieldWidth('%.2f',log_e{i}));
        qW = max(qW, findFieldWidth('%.2f',q{i}));
    end

    mW = findFieldWidth('%d',m);
    orderHeaderWidth = eW + qW + 1;

    fprintf('method: %s\nT: %d\n',methodName, T);

    % Print order headers
    fprintf(' %*s |',mW,'')
    for i = 1:length(orders)
        fprintf(' %-*s |', orderHeaderWidth, sprintf('Order %d', orders{i}));
    end
    fprintf('\n');


    % Print eq headers
    fprintf(' %*s |',mW,'m');
    for i = 1:length(orders)
        fprintf(' %*s %*s |', eW, 'e', qW, 'q');
    end
    fprintf('\n');


    % Print devider
    m_dev = repmat('-',1,mW);
    column_dev = repmat('-',1,orderHeaderWidth);
    fprintf('-%s-+',m_dev);
    for i = 1:length(orders)
        fprintf('-%s-+', column_dev);
    end
    fprintf('\n');



    % Print each row
    for i = 1:length(m)
        fprintf(' %*d |',mW,m(i));
        for j = 1:length(orders)
            if i == 1
                fprintf(' %*.2f %*s |', eW, log_e{j}(i), qW, '');
            else
                fprintf(' %*.2f %*.2f |', eW, log_e{j}(i), qW, q{j}(i-1));
            end
        end
        fprintf('\n');
    end

    fprintf('\n');

end