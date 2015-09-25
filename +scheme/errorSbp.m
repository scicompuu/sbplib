function e = errorSbp(discr, v1, v2)
    % If v1 and v2 are more complex types, something like grid functions... Then we may use .getVectorFrom here!
    H = discr.H;
    err = v2 - v1;
    e = sqrt(err'*H*err);
end