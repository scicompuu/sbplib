function alpha = matrixborrow(A,B)

    function y = neg(x)
        y = min(eig(A-x*B));
    end


    x0 = 1;
    while neg(x0) > -1
        x0 = 2*x0;
    end

    opt = optimset('Display','none');
    alpha = fsolve(@neg, x0,opt);
end
