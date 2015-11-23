function A = cell2sparse(C)

    if isempty(C)
        A = sparse([]);
        return
    end

    n = row_height(C);
    m = col_width(C);

    N = sum(n);
    M = sum(m);

    A = sparse(N,M);

    n_ind = [0 cumsum(n)];
    m_ind = [0 cumsum(m)];

    for i = 1:size(C,1)
        for j = 1:size(C,2)
            if ~has_matrix(C{i,j})
                continue
            end
            A(n_ind(i)+1:n_ind(i+1),m_ind(j)+1:m_ind(j+1)) = C{i,j};
        end
    end

end

function m = col_width(C)
    for j = 1:size(C,2)
        for i = 1:size(C,1)
            if ~has_matrix(C{i,j})
                continue
            end
            m(j) = size(C{i,j},2);
        end
    end
end

function n = row_height(C)
    for i = 1:size(C,1)
        for j = 1:size(C,2)
            if ~has_matrix(C{i,j})
                continue
            end
            n(i) = size(C{i,j},1);
        end
    end
end

function b = has_matrix(c)
    b = ~(isempty(c) || (numel(c)==1 && c == 0));
end