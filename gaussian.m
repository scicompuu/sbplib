function z = gaussian(x,x0,d)
    z = exp(-norm(x-x0).^2/d^2);
end