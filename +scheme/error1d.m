function e = error1d(schm, v1, v2)
    h = schm.h;
    e = sqrt(h*sum((v1-v2).^2));
end