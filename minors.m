function [minor, sub] = minors(A)
    [n, m] = size(A);

    if n ~= m
        error('A must be square');
    end

    sub = {};
    ks = {};

    ind = 1:n;
    for k = 1:n
        C = nchoosek(ind,k);
        for i = 1:size(C,1)
            ks{end + 1} = k;
            sub{end + 1} = A(C(i,:),C(i,:));
        end
    end

    for i = 1:length(sub)
        fprintf('%d:\n', ks{i});
        disp(sub{i})

        minor(i) = det(sub{i});
    end
end


% A is positive semidefinite if all minors are non-negative
% A is negative semidefinite if all odd minors are non-positive and all even minors are non-negative

