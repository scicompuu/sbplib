function print_issparse(A)
    warning('Deprecated! Use printExpr() instead!');
    b = issparse(A);
    if b
        s = 'true';
    else
        s = 'false';
    end
    fprintf('%8s is sparse: %s\n',inputname(1),s);
end