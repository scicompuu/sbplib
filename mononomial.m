function y = mononomial(x, k)
    if k < 0
        y = x*0;
        return
    end
    y = x.^k/factorial(k);
end

