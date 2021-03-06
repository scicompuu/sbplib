% Cell is a reimplementation of matlabs cell array with the benefit that it is subclassable
% It might be used for giving a typename to a cellarray to increase readability of the code.
classdef Cell
    properties
        data
    end

    methods
        function obj = Cell(data)
            default_arg('data', {});
            if ~iscell(data)
                error('Input argument to Cell must be a cell array');
            end

            obj.data = data;
        end

        function str = toString(obj)
            str = sprintf('%s%s', class(obj), toString(obj.data));
        end

        function s = size(A)
            s = size(A.data);
        end

        function b = isempty(A)
            b = prod(size(A)) == 0;
        end

        function l = length(A)
            l = length(A.data);
        end

        function ind = end(A,k,n)
            ind = builtin('end',A.data, k, n);
        end

        function B = transpose(A)
            b = A.data.';
            B = callConstructor(A, b);
        end

        function B = ctranspose(A)
            b = A.data';
            B = callConstructor(A, b);
        end

        function A = subsasgn(A, S, B)
            a = subsasgn(A.data, S, B);
            A = callConstructor(A, a);
        end

        function B = subsref(A, S)
            switch S(1).type
                case '()'
                    b = subsref(A.data, S(1));
                    B = callConstructor(A, b);
                    if length(S) > 1
                        B = subsref(B,S(2:end));
                    end
                case '{}'
                    B = subsref(A.data, S);
                case '.'
                    B = builtin('subsref',A, S);
                otherwise
                    error('unreachable');
            end
        end

        function C = horzcat(varargin)
            dataArray = cell(1, length(varargin));

            for i = 1:length(varargin)
                dataArray{i} = varargin{i}.data;
            end

            c = horzcat(dataArray{:});
            C = callConstructor(varargin{1}, c);
        end

        function C = vertcat(varargin)
            dataArray = cell(1, length(varargin));

            for i = 1:length(varargin)
                dataArray{i} = varargin{i}.data;
            end

            c = vertcat(dataArray{:});
            C = callConstructor(varargin{1}, c);
        end
    end
end
