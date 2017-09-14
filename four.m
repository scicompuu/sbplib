% four returns the fourier transform u_hat of the function u and the frequencies w
function [w, u_hat] = four(x, u)
    u_hat = fft(u);

    N = length(x);
    L = x(end) - x(1);

    k = shift_k(0:N-1);

    u_hat = fftshift(u_hat);

    dw = 2*pi/L;
    w = dw*k;
end

function k_shifted = shift_k(k)
    N = length(k);
    k_shifted = [-floor(N/2):-1, 0, 1:ceil(N/2)-1];
end
