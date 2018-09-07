function z = gaussian(x,x0,d)
    z = exp(-sum((x-x0).^2,2)/d^2);
end