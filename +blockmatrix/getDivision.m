function div = getDivision(bm)
    if ~blockmatrix.isBlockmatrix(bm)
        error('blockmatrix:getDivision:NotABlockmatrix', 'Input is not a blockmatrix');
    end

    if isempty(bm)
        div = {[],[]};
        return
    end

    div = {row_height(bm),col_width(bm)};
end


function m = col_width(C)
    m = zeros(1,size(C,2));
    for j = 1:size(C,2)
        for i = 1:size(C,1)
            if isNullMatrix(C{i,j})
                continue
            end
            m(j) = size(C{i,j},2);
        end
    end
end

function n = row_height(C)
    n = zeros(1,size(C,1));
    for i = 1:size(C,1)
        for j = 1:size(C,2)
            if isNullMatrix(C{i,j})
                continue
            end
            n(i) = size(C{i,j},1);
        end
    end
end

function b = isNullMatrix(A)
    [n, m] = size(A);
    b = n == 0 && m == 0;
end
