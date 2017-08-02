classdef Cell
    properties
        data
    end
    methods
        function obj = Cell(data)
            if ~iscell(data)
                class(data)
                error('Input argument to Cell must be a cell array');
            end

            obj.data = data;
        end

        % function display(A)
        %     n = size(A.data);

        %     sizeStr = join(cellfun(@num2str, num2cell(n), 'UniformOutput',false),'x');
        %     header = [sizeStr, 'Cell']

        %     disp()
        %     disp(A.data)
        %     % display(A.data)
        % end

        function disp(A)
            disp(A.data)
        end

        function A = subsasgn(A, S, B)
            a = subsasgn(A.data, S, B);
            A = callConstructor(A, a);
        end

        function B = subsref(A, S)
            switch S(1).type
                case '()'
                    b = subsref(A.data, S);
                    B = callConstructor(A, b);
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

        function vertcat(varargin)
            dataArray = cell(1, length(varargin));

            for i = 1:length(varargin)
                dataArray{i} = varargin{i}.data;
            end

            c = vertcat(dataArray{:});
            C = callConstructor(varargin{1}, c);
        end
    end
end
