% % Usage example:
% c = structCellArray({'a','b'}, {
%     1, 2;
%     3, 4;
% });

function c = structCellArray(fields, values)
    assert(length(fields) == size(values, 2), 'Number of fields and number of colums of ''values'' must be equal');
    c = cell(1, size(values, 1));

    for i = 1:size(values, 1)
        for j = 1:length(fields)
            c{i}.(fields{j}) = values{i,j};
        end
    end
end
