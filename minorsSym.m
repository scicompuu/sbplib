function ineq = minorsSym(B)
    B

    [n, m] = size(B);

    if n ~= m
        error('B must be square');
    end

    subB = {};
    ks = {};

    ind = 1:n;
    for k = 1:n
        C = nchoosek(ind,k);
        for i = 1:size(C,1)
            ks{end + 1} = k;
            subB{end + 1} = B(C(i,:),C(i,:));
        end
    end



    for i = 1:length(subB)
        fprintf('%d:\n', ks{i});
        disp(subB{i})

        minor{i} = det(subB{i});
    end

    for i = 1:length(subB)
        fprintf('%d:\n', ks{i});
        disp(minor{i});

        if mod(ks{i},2) == 0
            ineq{i} = minor{i} >= 0;
        else
            ineq{i} = minor{i} <= 0;
        end
    end

    ineqsys = true;
    for i = 1:length(ineq)
        ineqsys = ineqsys & ineq{i};
    end

    ineq = simplify(ineqsys)
end


% B is positive semidefinite if all minors are non-negative
% B is negative semidefinite if all odd minors are non-positive and all even minors are non-negative

