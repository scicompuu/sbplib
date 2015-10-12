function printSize(A)
    warning('Deprecated! Use printExpr() instead!');
    s = size(A);
    fprintf('%8s has size: [%d, %d]\n',inputname(1),s(1),s(2));
end