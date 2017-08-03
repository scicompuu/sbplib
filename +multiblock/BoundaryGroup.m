% BoundaryGroup defines a boundary grouping in a multiblock grid.
classdef BoundaryGroup < Cell
    methods
        function obj = BoundaryGroup(data)
            obj = obj@Cell(data);
        end

        % function display(obj, name)

        %     disp(' ')
        %     disp([name, ' ='])
        %     disp(' ')

        %     if length(obj.names) == 1
        %         fprintf('    {}\n\n')
        %         return
        %     end

        %     fprintf('    {')

        %     fprintf('%d:%s', obj.blockIDs(1), obj.names{1})
        %     for i = 2:length(obj.names)
        %         fprintf(', %d:%s', obj.blockIDs(i), obj.names{i});
        %     end

        %     fprintf('}\n\n')
        % end
    end
end