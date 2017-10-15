function u = ifour(u_hat)
    u_hat = ifftshift(u_hat);
    u = ifft(u_hat);
end
