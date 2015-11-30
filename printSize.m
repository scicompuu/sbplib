function printSize(A)
    result = size(A);
    fprintf('size(%s) => %s\n', inputname(1), toString(result));
end