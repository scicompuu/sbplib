function e = errorRelative(~,v1,v2)
    e = sqrt(sum((v1-v2).^2)/sum(v2.^2));
end