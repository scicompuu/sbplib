function bm = fromMatrix(A, div)
    d1 = div{1};
    d2 = div{2};
    [n, m] = size(A);
    if n ~= sum(d1) || m ~= sum(d2)
        error('blockmatrix:fromMatrix:NonMatchingDim','The dimensions in div does not sum to the dimensions in A.');
    end

    bm = cell(length(d1), length(d2));
    I = 1;
    for i = 1:length(d1)
        J = 1;
        for j = 1:length(d2)
            Asub = A(I:(I + d1(i)-1), J:(J + d2(j)-1));
            if nnz(Asub) == 0
                bm{i,j} = [];
            else
                bm{i,j} = Asub;
            end
            J = J + d2(j);
        end
        I = I + d1(i);
    end
end