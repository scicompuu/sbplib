% Calculates the maximum field width needed width fprintf for a given format fmt for the values in A
%     A = [1.13232 10.233 1.1];
%     width = findFieldWidth('%.2f',A);
% Gives wdith = 5
function width = findFieldWidth(fmt, A)
    width = 0;
    for i = 1:numel(A)
        elem = A(i);
        if iscell(elem)
            elem = elem{1};
        end
        str = sprintf(fmt,elem);
        width = max(width,length(str));
    end
end
