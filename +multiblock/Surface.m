classdef Surface < handle
    properties
        grid
        surfs

        ZData
        CData

    end

    methods
        function obj = Surface(g, gf)
            obj.grid = g;

            coords = obj.grid.points();
            X = obj.grid.funcToPlotMatrices(coords(:,1));
            Y = obj.grid.funcToPlotMatrices(coords(:,2));

            V = obj.grid.funcToPlotMatrices(gf);


            holdState = ishold();
            hold on

            surfs = cell(1, obj.grid.nBlocks);
            for i = 1:obj.grid.nBlocks
                surfs{i} = surf(X{i}, Y{i}, V{i});
            end

            if holdState == false
                hold off
            end

            obj.surfs = [surfs{:}];

            obj.ZData = gf;
            obj.CData = gf;
        end

        function set(obj, propertyName, propertyValue)
            set(obj.surfs, propertyName, propertyValue);
        end

        function obj = set.ZData(obj, gf)
            obj.ZData = gf;

            V = obj.grid.funcToPlotMatrices(gf);
            for i = 1:obj.grid.nBlocks
                obj.surfs(i).ZData = V{i};
            end
        end

        function obj = set.CData(obj, gf)
            obj.CData = gf;

            V = obj.grid.funcToPlotMatrices(gf);
            for i = 1:obj.grid.nBlocks
                obj.surfs(i).CData = V{i};
            end
        end
    end
end
