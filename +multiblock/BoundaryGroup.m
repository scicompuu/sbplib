% BoundaryGroup defines a boundary grouping in a multiblock grid.
classdef BoundaryGroup
    properties
        blockIDs
        names
    end

    methods
        function obj = BoundaryGroup(varargin)
            % Input arguemnts are arbitrary number or 1x2 cell arrays
            % representing each boundary in the group.
            % The 1st element of the cell array is an integer defining which grid it belongs to.
            % The 2nd element of the cell array is the name of the boundary within the block.
            %
            % Ex:
            %   bg = multiblock.BoundaryGroup({1,'n'},{1,'s'},{2,'s'})


            obj.blockIDs = [];
            obj.names = {};
            for i = 1:length(varargin)
                obj.blockIDs(i) = varargin{i}{1};
                obj.names{i} = varargin{i}{2};
            end
        end

        function display(obj, name)

            disp(' ')
            disp([name, ' ='])
            disp(' ')

            if length(obj.names) == 1
                fprintf('    {}\n\n')
                return
            end

            fprintf('    {')

            fprintf('%d:%s', obj.blockIDs(1), obj.names{1})
            for i = 2:length(obj.names)
                fprintf(', %d:%s', obj.blockIDs(i), obj.names{i});
            end

            fprintf('}\n\n')
        end
    end
end