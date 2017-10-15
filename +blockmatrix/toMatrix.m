function A = toMatrix(bm)
    if ~blockmatrix.isBlockmatrix(bm)
        error('blockmatrix:toMatrix:NotABlockmatrix', 'Input is not a blockmatrix');
    end

    div = blockmatrix.getDivision(bm);
    n = div{1};
    m = div{2};

    N = sum(n);
    M = sum(m);

    A = sparse(N,M);

    n_ind = [0 cumsum(n)];
    m_ind = [0 cumsum(m)];

    for i = 1:size(bm,1)
        for j = 1:size(bm,2)
            if(isempty(bm{i,j}))
                bm{i,j} = sparse(n(i),m(j));
            end
        end
    end
    A = cell2mat(bm);
end
