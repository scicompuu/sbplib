% TODO
%   Documentation
%   Should the storage be done in some other way?
%   Cell array as keys?
%   Some possibility to load and save to a matfile?
%       May be load and save from outside is enough.
%   Get and set functions called by subsref and subsasgn?
%   Unit tests.


classdef Dictionary
    properties
        store
    end

    methods
        function obj = Dictionary()
            obj.store = struct();
        end

        function s = getStore(obj)
            s = obj.store;
        end

        function display(obj)
            if length(fieldnames(obj.store)) == 0
                fprintf('%s is an empty Dictionary\n',inputname(1));
                return
            end

            lineformat = [inputname(1) '(%s) = %s\n'];

            display_impl(obj.store,'');

            function display_impl(s, path)
                if ~isstruct(s)
                    fprintf(lineformat,path(3:end), value2str(s));
                    % fprintf(') = %s\n', value2str(s));
                    % fprintf('%s(', objName);
                    return
                end

                fn = fieldnames(s);

                for i = 1:length(fn)
                    display_impl(s.(fn{i}), [path ', ' fn{i}(2:end)]);
                end

            end

            function str = value2str(val)
                if isnumeric(val) || ischar(val)
                    str = mat2str(val);
                else
                    str = class(val);
                end
            end
        end

        % B = obj(i)
        function B = subsref(obj,S)
            switch S.type
                case '()'
                    Sf = obj.subs2dotSubs(S);
                    try
                        B = subsref(obj.store,Sf);
                    catch ME
                        if strcmp(ME.identifier,'MATLAB:nonExistentField')
                            error('Reference to non-existent entry %s',toString(S.subs));
                        else
                            throw(ME);
                        end
                    end
                otherwise
                    B = builtin('subsref', obj, S);
                    % error('Unsupported indexing operator: %s',S.type);
            end
        end

         % A(i) = B
        function obj = subsasgn(obj,S,B);
            switch S.type
                case '()'
                    Sf = obj.subs2dotSubs(S);
                    obj.store = subsasgn(obj.store,Sf,B);
                otherwise
                    error('Unsupported indexing operator: %s',S.type);
            end
        end

        function Sf = subs2dotSubs(obj,S)
            for i = 1:length(S.subs)
                Sf(i).type = '.';
                Sf(i).subs = obj.getFieldname(S.subs{i});
            end
        end

        % Should probably use mat2str with some kind of normalization to make all variables valied fieldname
        %  and make it possible to recover the value
        function fName = getFieldname(obj, val)
            if isnumeric(val)
                valStr = num2str(val);
            elseif ischar(val)
                valStr = val;
            else
                error('Dont know what to do with val!');
            end
            fName = sprintf('f%s',valStr);
        end
    end
end