% Creates a block matrix according to the division with zeros everywhere.
function bm = zero(div)
    n = div{1};
    m = div{2};

    N = length(n);
    M = length(m);

    bm = cell(N,M);

    for i = 1:N
        for j = 1:M
            bm{i,j} = sparse(n(i),m(j));
        end
    end
end
