function A = diags(B,d,m,n)
    assert(size(B,1) == m);

    A = repmat(B(:,1)*0, [1, n]);

    for i = 1:size(B,2)
        A(:,d(i)+ (1:m)) = A(:,d(i)+ (1:m)) + diag(B(:,i));
    end
end