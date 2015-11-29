function I = fzero_vec(f,lim)
% Wrapper around the built-in function fzero that
% handles multiple functions (vector-valued f).

    fval = f(lim(1));
    I = zeros(size(fval));
    
    for i = 1:length(fval)
        e = zeros(size(fval));
        e(i) = 1;
        I(i) = fzero(@(t) sum(e.*f(t)),lim);
    end
    
end