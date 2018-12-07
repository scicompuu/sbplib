classdef Map < handle
    properties
        map
    end

    % can we support multi map using varargin?
    % probably a bad idea. For example it complicates keys();

    methods
        function obj = Map()
            obj.map = containers.Map();
        end

        function set(obj, k, v)
            keyByteStream = getByteStreamFromArray(k);

            obj.map(char(keyByteStream)) = v;
        end

        function v = get(obj, k)
            keyByteStream = getByteStreamFromArray(k);

            v = obj.map(char(keyByteStream));
        end

        function b = isKey(obj, k)
            keyByteStream = getByteStreamFromArray(k);
            b = obj.map.isKey(char(keyByteStream));
        end

        function c = keys(obj)
            keyByteStreams = obj.map.keys;

            n = length(keyByteStreams);

            c = cell(1, n);
            for i = 1:n
                c{i} = getArrayFromByteStream(uint8(keyByteStreams{i}));
            end
        end

        function l = length(obj)
            l = obj.map.length;
        end

        function remove(obj, k)
            keyByteStream = getByteStreamFromArray(k);
            obj.map.remove(char(keyByteStream));
        end

        function s = size(obj)
            s = obj.map.size;
        end

        function c = values(obj)
            c = obj.map.values;
        end

        function v = subsref(obj, S)
            switch S(1).type
                case '()'
                    if length(S.subs) > 1
                        error('sbplib:Map:multipleKeys', 'Multiple dimensions are not supported. Use a cell array as a key instead.');
                    end
                    k = S.subs{1};
                    try
                        v = get(obj, k);
                    catch ME
                        if strcmp(ME.identifier,'MATLAB:Containers:Map:NoKey')
                            error('Reference to non-existent entry %s',toString(S.subs));
                        else
                            throw(ME);
                        end
                    end
                otherwise
                    try
                        v = builtin('subsref', obj, S);
                    catch e
                        error('You can''t use dot notation for this because Matlab(TM). What is this piece of shit software anyway?')
                    end
            end
        end

        function obj = subsasgn(obj, S, v);
            switch S(1).type
                case '()'
                    if length(S.subs) > 1
                        error('sbplib:Map:multipleKeys', 'Multiple dimensions are not supported. Use a cell array as a key instead.');
                    end
                    k = S.subs{1};
                    set(obj, k, v);
                otherwise
                    error('You can''t use dot notation because Matlab(TM). What is this piece of shit software anyway?')
            end
        end
    end
end
