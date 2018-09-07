% BoundaryGroup defines a boundary grouping in a multiblock grid.
% It workds like a cell array and collects boundary identifiers
% Within the multiblock package a BoundaryGroup is a valid boundary identifier as well.
classdef BoundaryGroup < Cell
    methods
        function obj = BoundaryGroup(data)
            obj = obj@Cell(data);
        end

        function display(obj, name)

            disp(' ')
            disp([name, ' ='])
            disp(' ')

            fprintf('    BoundaryGroup%s\n\n', toString(obj.data));
        end
    end
end
