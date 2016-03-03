% sparse2cell breaks a sparse matrix up into a cell matrix of sparse matrices.
% any zero submatrix creates a empty cell in the cell matrix.
%  A -- NxM sparse matrix
%  d1, d2 -- vectors of sub matrix sizes for each dimensions. Must have sum(di) == Ni.
% Example:
%   C = sparse2cell(A,[5 10], [10 5])
function C = sparse2cell(A, d1, d2)
    [n, m] = size(A);
    if n ~= sum(d1) || m ~= sum(d2)
        error('sparse2cell:NonMatchingDim','The elements of d1 and d2 must sum to N and M.');
    end

    C = cell(length(d1), length(d2));
    I = 1;
    for i = 1:length(d1)
        J = 1;
        for j = 1:length(d2)
            Asub = A(I:(I + d1(i)-1), J:(J + d2(j)-1));
            if nnz(Asub) == 0
                C{i,j} = [];
            else
                C{i,j} = Asub;
            end
            J = J + d2(j);
        end
        I = I + d1(i);
    end
end
