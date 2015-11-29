function I = integral_vec(f,a,b)
% Wrapper around the built-in function integral that
% handles multiple limits.

    Na = length(a);
    Nb = length(b);
    assert(Na == 1 || Nb == 1 || Na==Nb,...
        'a and b must have same length, unless one is a scalar.');
    
    if(Na>Nb);
        I = zeros(size(a));
        for i = 1:Na
            I(i) = integral(f,a(i),b);
        end
    elseif(Nb>Na)
        I = zeros(size(b));
        for i = 1:Nb
            I(i) = integral(f,a,b(i));
        end
    else
        I = zeros(size(b));
        for i = 1:Nb
            I(i) = integral(f,a(i),b(i));
        end
    end
end