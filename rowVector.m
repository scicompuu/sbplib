function r = rowVector(v)
    if  size(v,1) == 1
        r = v
    elseif size(v,2) == 1
        r = v';
    else
        error('v is not a matrix');
    end
end