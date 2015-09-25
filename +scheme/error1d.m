function e = error1d(discr, v1, v2)
    h = discr.h;
    e = sqrt(h*sum((v1-v2).^2));
end