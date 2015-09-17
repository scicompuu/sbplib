% Replaces string old on the console with new.
% new may be a format string like the one accepted by fprintf
% old has to be the last string printed.
function s = replace_string(old,new,varargin)
    reverseStr = repmat(sprintf('\b'), 1, length(old));
    blankStr   = repmat(sprintf(' ') , 1, length(old));
    fprintf([reverseStr,blankStr,reverseStr, new],varargin{:});
    s = sprintf(new, varargin{:});
end