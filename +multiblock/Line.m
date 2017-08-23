classdef Line < handle
    properties
        grid
        lines

        YData
    end

    methods
        function obj = Line(g, gf)
            assertType(g, 'multiblock.Grid')
            obj.grid = g;

            X = obj.grid.splitFunc(obj.grid.points());
            Y = obj.grid.splitFunc(gf);

            holdState = ishold();
            hold on

            lines = {1, obj.grid.nBlocks};
            for i = 1:obj.grid.nBlocks
                lines{i} = plot(X{i}, Y{i});
            end

            if holdState == false
                hold off
            end

            obj.lines = [lines{:}];

            obj.YData = gf;
        end

        function set(obj, propertyName, propertyValue)
            set(obj.lines, propertyName, propertyValue);
        end

        function obj = set.YData(obj, gf)
            obj.YData = gf;

            Y = obj.grid.funcToPlotMatrices(gf);
            for i = 1:obj.grid.nBlocks
                obj.lines(i).YData = Y{i};
            end
        end
    end
end
