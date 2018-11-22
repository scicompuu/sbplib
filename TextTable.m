classdef TextTable < handle
    properties
        data
        fmtArray
        vertDiv
        horzDiv
    end

    methods
        function obj = TextTable(data, vertDiv, horzDiv)
            default_arg('vertDiv', []);
            default_arg('horzDiv', []);

            obj.data = data;
            obj.vertDiv = vertDiv;
            obj.horzDiv = horzDiv;

            obj.fmtArray = cell(size(data));
            obj.formatAll('%s');

        end

        function n = nRows(obj)
            n = size(obj.data, 1);
        end

        function n = nCols(obj)
            n = size(obj.data, 2);
        end

        function print(obj)
            disp(obj.toString());
        end

        function formatAll(obj, fmt)
            obj.fmtArray = cell(size(obj.data));
            obj.fmtArray(:,:) = {fmt};
        end

        function formatCell(obj, i, j, fmt)
            obj.fmtArray{i,j} = fmt;
        end

        function formatGrid(obj, I, J, fmt)
            for i = I
                for j = J
                    obj.fmtArray{i,j} = fmt;
                end
            end
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

            N = size(strArray, 2);

            % First horzDiv
            if isDiv(0, obj.horzDiv, N);
                str = [str, obj.getHorzDiv(widths)];
            end

            for i = 1:obj.nRows
                str = [str, TextTable.rowToString(strArray(i,:), widths, obj.vertDiv)];

                % Interior horzDiv
                if isDiv(i, obj.horzDiv, N)
                    str = [str, obj.getHorzDiv(widths)];
                end
            end
        end

        function str = getHorzDiv(obj, widths)
            str = TextTable.rowToString(cell(1,obj.nCols), widths, obj.vertDiv);
            str(find(' ' == str)) = '-';
            str(find('|' == str)) = '+';
        end

        function strArray = getStringArray(obj)
            assert(all(size(obj.data) == size(obj.fmtArray)), 'Sizes of format matrix and data matrix do not match.')
            strArray = cell(size(obj.data));

            for i = 1:obj.nRows
                for j = 1:obj.nCols
                    strArray{i,j} = sprintf(obj.fmtArray{i,j}, obj.data{i,j});
                end
            end
        end
    end

    methods (Static)
        function str = rowToString(strs, widths, vertDiv)
            N = length(strs);

            % First vertDiv
            if isDiv(0, vertDiv, N)
                prefix = '| ';
            else
                prefix = ' ';
            end

            % Pad strings
            for i = 1:N
                strs{i} = sprintf('%*s', widths(i), strs{i});
            end

            % Column delimiters
            delims = cell(1,N-1);
            for i = 1:length(delims)
                if isDiv(i, vertDiv, N);
                    delims{i} = '| ';
                else
                    delims{i} = ' ';
                end
            end

            if isDiv(N, vertDiv, N);
                suffix = '|';
            else
                suffix = '';
            end

            str = [prefix, strjoin(strs, delims), suffix, sprintf('\n')];
        end
    end
end

function b = isDiv(i, div, N)
    b = ismember(i, div) || ismember(i, N+div+1);
end