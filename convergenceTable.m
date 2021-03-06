function convergenceTable(caption, orders, m, e, q, tableType)
    default_arg('tableType','plaintext')

    switch tableType
        case {'plaintext','text','plain'}
            plainTextTable(caption, orders, m, e, q);
        case {'tex', 'latex'}
            latexTable(caption, orders, m, e, q);
    end
end

function plainTextTable(caption, orders, m, e, q)


    eW = 0;
    qW = 0;
    for i = 1:length(orders)
        log_e{i} = log10(e{i});
        eW = max(eW, findFieldWidth('%.2f',log_e{i}));
        qW = max(qW, findFieldWidth('%.2f',q{i}));
    end

    mW = findFieldWidth('%d',m);
    orderHeaderWidth = eW + qW + 1;

    fprintf('%s\n',caption);

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

function latexTable(caption, orders, m, e, q)

    nOrders = length(orders);

    header = {
        '\begin{table}[H]'
        '\centering'
        ['\begin{tabular}{c' repmat('|cc',1,nOrders) '} &']
        orderheaders(orders)
        '\hline'
        ['$N$'   repmat('& $log_{10}(l_2)$ & $q$',1,nOrders) ' \\']
        '\hline'
    };

    footer = {
        '\end{tabular}'
        ['\caption{' caption '}']
        '\label{table:LABEL}'
        '\end{table}'
    };

    data = cell(1,length(m));
    data{1} = num2str(m(1));
    for j = 1:nOrders
        data{1} = [data{1} ' & ' sprintf('%8.2f',log10(e{j}(1))) ' &         ' ];
    end
    data{1} = [data{1} '\\'];

    for i = 2:length(m)
        data{i} = [data{i} num2str(m(i))  ];
        for j = 1:nOrders
            data{i} = [data{i} ' & ' sprintf('%8.2f',log10(e{j}(i)))  ' & '  sprintf('%8.2f',(q{j}(i-1))) ];
        end
        data{i} = [data{i} '\\'];
    end

    nlc = sprintf('\n');

    header = strjoin(header', nlc);
    data = strjoin(data, nlc);
    footer = strjoin(footer', nlc);

    table = strjoin({header, data, footer}, nlc);
    fprintf('%s\n', table);
end



function s = orderheaders(orders)
    s= sprintf('\\multicolumn{2}{|c}{%dth order}',orders{1});
    nOrders = length(orders);
    for i = 2:nOrders
        s = [s sprintf('& \\multicolumn{2}{|c}{%dth order}',orders{i})];
    end
    s = [s ' \\'];
end
