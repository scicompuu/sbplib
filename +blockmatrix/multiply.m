function C = multiply(A, B)
    [n_a, m_a] = size(A);
    [n_b, m_b] = size(B);

    assert(m_a == n_b, 'Dimensions do not agree.')

    C = cell(n_a, m_b);

    for i = 1:n_a
        for j = 1:m_b
            C{i,j} = sparse(size(A{i,1},1), size(B{1,j},2));
            for k = 1:n_b
                C{i,j} = C{i,j} + A{i,k}*B{k,j};
            end
        end
    end
end
