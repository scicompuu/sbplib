function b = isBlockmatrix(bm)
    if ~iscell(bm)
        b = false;
        return
    end

    % Make sure all blocks are numerical matrices
    for i = 1:length(bm)
        if ~isnumeric(bm{i})
            b = false;
            return
        end
    end

    [N,M] = size(bm);
    % Make sure column dimensions agree
    for i = 1:N
        d = [];
        for j = 1:M
            d_ij = size(bm{i,j},1);
            if d_ij == 0
                continue
            end

            if isempty(d)
                d = d_ij;
                continue
            end

            if d ~= d_ij
                b = false;
                return
            end
        end
    end

    % Make sure row dimensions agree
    for j = 1:M
        d = [];
        for i = 1:N
            d_ij = size(bm{i,j},2);
            if d_ij == 0
                continue
            end

            if isempty(d)
                d = d_ij;
                continue
            end

            if d ~= d_ij
                b = false;
                return
            end
        end
    end

    b = true;
end
