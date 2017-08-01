classdef TextTable < handle
    properties
        data
        fmtArray
        vertDiv
        horzDiv

        nCols
        nRows
    end

    methods
        function obj = TextTable(data, vertDiv, horzDiv);
            default_arg('vertDiv', []);
            default_arg('horzDiv', []);


            obj.data = data;
            obj.vertDiv = vertDiv;
            obj.horzDiv = horzDiv;

            [obj.nRows, obj.nCols] = size(data);
            obj.fmtArray = cell(size(data));
            obj.formatAll('%s');

        end

        function formatAll(obj, fmt)
            obj.fmtArray(:,:) = {fmt};
        end

        function formatCell(obj, i, j, fmt)
            obj.fmtArray{i,j} = fmt;
        end

        function formatRow(obj, i, fmt)
            obj.fmtArray(i,:) = {fmt};
        end

        function formatColumn(obj, j, fmt)
            obj.fmtArray(:,j) = {fmt};
        end

        function widths = getWidths(obj)
            strArray = obj.getStringArray();

            widths = zeros(1,obj.nCols);
            for j = 1:obj.nCols
                for i = 1:obj.nRows
                    widths(j) = max(widths(j), length(strArray{i,j}));
                end
            end
        end

        function str = toString(obj)
            strArray = obj.getStringArray();
            widths = obj.getWidths();

            str = '';

            % First horzDiv
            if ismember(0, obj.horzDiv)
                str = [str, obj.getHorzDiv(widths)];
            end

            for i = 1:obj.nRows
                str = [str, TextTable.rowToString(strArray(i,:), widths, obj.vertDiv, obj.horzDiv)];

                % Interior horzDiv
                if ismember(i, obj.horzDiv)
                    str = [str, obj.getHorzDiv(widths)];
                end
            end
        end

        function str = getHorzDiv(obj, widths)
            str = TextTable.rowToString(cell(1,obj.nCols), widths, obj.vertDiv, obj.horzDiv);
            str(find(' ' == str)) = '-';
            str(find('|' == str)) = '+';
        end

        function strArray = getStringArray(obj)
            strArray = cell(size(obj.data));

            for i = 1:obj.nRows
                for j = 1:obj.nCols
                    strArray{i,j} = sprintf(obj.fmtArray{i,j}, obj.data{i,j});
                end
            end
        end
    end

    methods (Static)
        function str = rowToString(strs, widths, vertDiv, horzDiv)
            % First vertDiv
            if ismember(0, vertDiv)
                str = '| ';
            else
                str = ' ';
            end

            % Interior cols
            for j = 1:length(strs) - 1
                str = [str, sprintf('%*s ', widths(j), strs{j})];

                % Interior vertDiv
                if ismember(j, vertDiv)
                    str = [str, '| '];
                end
            end

            % Last col
            str = [str, sprintf('%*s ', widths(end), strs{end})];

            if ismember(length(strs), vertDiv)
                str = [str, '|'];
            end

            str = [str, sprintf('\n')];
        end
    end
end