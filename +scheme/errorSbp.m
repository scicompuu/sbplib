function e = errorSbp(schm, v1, v2)
    H = schm.H;
    err = v2 - v1;
    e = sqrt(err'*H*err);
end