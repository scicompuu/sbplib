% Copy the struct src to dest with default values from default
%   dest = copyWithDefault(src, default)
function dest = copyWithDefault(src, default)
    % src does not have a value => use default
    if isempty(src)
        dest = default;
        return
    end

    % src has a value and is not a struct => use src
    % src has a value and default is not a struct => use src
    if ~isstruct(src) || ~isstruct(default)
        dest = src;
        return
    end


    % src has a value and is a struct => add all default fields
    dest = src;

    fn = fieldnames(default);
    for i = 1:length(fn)
        if isfield(src, fn{i})
            srcField = src.(fn{i});
        else
            srcField = [];
        end

        dest.(fn{i}) = copyWithDefault(srcField, default.(fn{i}));
    end
end
