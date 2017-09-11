function b = isDivision(div)
    % Make sure it is a cellarray
    if ~iscell(div)
        b = false;
        return
    end

    % Make sure it has the right shape
    if numel(div) ~= 2
        b = false;
        return
    end

    if ~isDivisionVector(div{1}) || ~isDivisionVector(div{2})
        b = false;
        return
    end

    b = true;
end

function b = isDivisionVector(v)
    if isempty(v)
        b = true;
        return
    end

    if any(v < 0)
        b = false;
        return
    end

    b = true;
end
