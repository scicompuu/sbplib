function e = error2d(discr, v1, v2)
    % If v1 and v2 are more complex types, something like grid functions... Then we may use .getVectorFrom here!
    h = discr.h;
    e = sqrt(h.^2*sum((v1-v2).^2));
end