function e = error2d(schm, v1, v2)
    hu = schm.h(1);
    hv = schm.h(2);
    e = sqrt(hu*hv*sum((v1-v2).^2));
end