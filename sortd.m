% Sorts x with respect to the function d
% d is a function that accepts an element of x and returns a scalar.
function y = sortd(d,x)
    v = [];
    for i = 1:length(x)
        v(i) = d(x(i));
    end

    [~,I] = sort(v);
    y = x(I);
end