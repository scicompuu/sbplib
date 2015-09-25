function q = convergence(e,h)
    for i = 1:length(e)-1
        q(i) = log(e(i)/e(i+1))/log(h(i+1)/h(i));
    end
end