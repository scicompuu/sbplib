function A = spdiagVariable(a,i)
    default_arg('i',0);

    if isrow(a)
        a = a';
    end

    n = length(a)+abs(i);

    if i > 0
    	a = [sparse(i,1); a];
    elseif i < 0
    	a = [a; sparse(abs(i),1)];
    end

    A = spdiags(a,i,n,n);
end