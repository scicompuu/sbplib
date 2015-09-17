% I = find_elements(a,b)
% Finds the index of elements a in b.
% a and b have to be in the same order.
function I = find_elements(a,b)
    I = [];

    j = 1;
    for i = 1:length(a)
        while b(j) ~= a(i)
            j = j + 1;
        end
        I(end+1) = j;
        j = j+1;
    end

    assert(length(I) == length(a),'Expected %d but got %d elements',length(a),length(I))
end