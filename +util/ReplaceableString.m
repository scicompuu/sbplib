classdef ReplaceableString < handle
    properties
        n
        fmt
        param
    end

    methods
        function obj = ReplaceableString(fmt, varargin)
            default_arg('fmt', '');
            obj.n = 0;

            obj.setFormat(fmt);
            obj.param = varargin;
        end

        function setFormat(obj, fmt)
            obj.fmt = fmt;
        end

        function update(obj, fmt, varargin)
            obj.setFormat(fmt);
            obj.param = varargin;

            obj.display();
        end

        function updateParam(obj, varargin)
            obj.param = varargin;
            obj.display();
        end

        function display(obj)
            reverseStr = repmat(sprintf('\b'), 1, obj.n);
            newStr = padStr(sprintf(obj.fmt, obj.param{:}),obj.n);
            fprintf([reverseStr, newStr]);

            obj.n = length(newStr);
        end
    end

end

function b = padStr(a, n)
    b = sprintf('%-*s', n, a);
end
