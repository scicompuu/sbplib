classdef Contour < handle
    properties
        grid
        contours
        nContours

        ZData
        CData

    end

    methods
        function obj = Contour(g, gf, nContours)
            obj.grid = g;
            obj.nContours = nContours;

            coords = obj.grid.points();
            X = obj.grid.funcToPlotMatrices(coords(:,1));
            Y = obj.grid.funcToPlotMatrices(coords(:,2));

            V = obj.grid.funcToPlotMatrices(gf);


            holdState = ishold();
            hold on

            contours = {1, obj.grid.nBlocks};
            for i = 1:obj.grid.nBlocks
                [~, contours{i}] = contour(X{i}, Y{i}, V{i},obj.nContours);
                contours{i}.LevelList = contours{1}.LevelList;
            end

            if holdState == false
                hold off
            end

            obj.contours = [contours{:}];

            obj.ZData = gf;
            obj.CData = gf;
        end

        function set(obj, propertyName, propertyValue)
            set(obj.contours, propertyName, propertyValue);
        end

        function obj = set.ZData(obj, gf)
            obj.ZData = gf;

            V = obj.grid.funcToPlotMatrices(gf);
            for i = 1:obj.grid.nBlocks
                obj.contours(i).ZData = V{i};
            end
        end

        function obj = set.CData(obj, gf)
            obj.CData = gf;

            V = obj.grid.funcToPlotMatrices(gf);
            for i = 1:obj.grid.nBlocks
                obj.contours(i).CData = V{i};
            end
        end
    end
end
