function y = rickerWavelet(x, x0, A)
    y = (1-2*pi^2*A^2*(x-x0).^2).*exp(-pi^2*A^2*(x-x0).^2);
end
