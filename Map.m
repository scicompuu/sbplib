classdef Map < handle
    properties
        map
    end

    % can we support multi map using varargin?

    methods
        function obj = Map()
            obj.map = containers.Map()
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
            b = obj.map.isKey(keyByteStream);
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
            switch S.type
                case '()'
                    k = S.subs;
                    try
                        v = obj.get(k);
                    catch ME
                        if strcmp(ME.identifier,'MATLAB:Containers:Map:NoKey')
                            error('Reference to non-existent entry %s',toString(S.subs));
                        else
                            throw(ME);
                        end
                    end
                otherwise
                    v = builtin('subsref', obj, S);
            end
        end

        function obj = subsasgn(obj, S, v);
            switch S.type
                case '()'
                    k = S.subs;
                    obj.set(k, v);
                otherwise
                    error('Unsupported indexing operator: %s',S.type);
            end
        end
    end
end
