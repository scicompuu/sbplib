classdef SolutionFile < handle
    properties
        filename
        matfile
        keys % Cell array of keys. Each key is a structure
        % entries
    end

    methods
        function obj = SolutionFile(filename)

            obj.filename = filename;

            is_new_file = ~exist(filename,'file');

            % obj.matfile = matfile(filename,'Writable',true);
            fprintf('MATLAB SUCKS!!!!\n')

            if is_new_file
                obj.matfile.keys = {};
                obj.matfile.entries = {};
            else
                matObj = matfile(filename,'Writable',true);
                obj.matfile.keys = matObj.keys;
                obj.matfile.entries = matObj.entries;
            end

            obj.keys = obj.matfile.keys;

        end

        function stupidSave(obj)
            matObj = matfile(obj.filename,'Writable',true);

            keys = obj.matfile.keys;
            entries = obj.matfile.entries;

            delete(obj.filename);

            matObj = matfile(obj.filename,'Writable',true);
            matObj.keys = keys;
            matObj.entries = entries;
        end

        function list(obj, show_syntax)
            default_arg('show_syntax',false);
            for i = 1:length(obj.keys)
                fprintf('[%d]: %s', i, struct2string(obj.keys{i}));

                if show_syntax
                    fprintf('\t%s',struct2syntax(obj.keys{i}));
                end

                fprintf('\n');
            end
        end

        % Returns the index for a key. Returns 0 if the key doesn't exist
        function I = getIndex(obj, key)
            I = 0;

            for i = 1:length(obj.keys)
                if isequal(key,obj.keys{i});
                    I = i;
                    return
                end
            end
        end

        function b = isKey(obj, key)
            I = obj.getIndex(key);

            if I == 0
                b = false;
            else
                b = true;
            end
        end

        % Gets entries where key match exactly.
        function e = get(obj, key)
            if ~obj.isKey(key);
                error('No such key: %s', struct2string(key));
            end

            I = obj.getIndex(key);
            e = obj.getEntryByIndex(I); % unpack the cell array
        end


        % Handles indexing weirdness of matfile class
        function e = getEntryByIndex(obj, I)
            e = obj.matfile.entries(1,I);
            e = e{1};
        end

        function e = deleteByIndex(obj, I)
            obj.keys(I) = [];
            obj.matfile.keys = obj.keys;

            entries = obj.matfile.entries;
            entries(I) = [];
            obj.matfile.entries = entries;
        end

        % Handles indexing weirdness of matfile class
        function setEntryByIndex(obj,I, entry, force_flag)
            default_arg('force_flag',false);
            if ~force_flag && ( I < 1 || I > length(obj.keys))
                error('I is out of range. I = %d',I);
            end
            obj.matfile.entries(1,I) = {entry};
        end

        function store(obj, key, entry)
            if obj.isKey(key);
                I = obj.getIndex(key);
            else
                I = length(obj.keys) + 1;
            end

            obj.keys{I} = key;
            obj.matfile.keys = obj.keys;
            obj.setEntryByIndex(I,entry,true);
        end

        function delete(obj, key)
            if ~obj.isKey(key);
                error('No such key: %s', struct2string(key));
            end
            I = obj.getIndex(key);
            obj.deleteByIndex(I);
        end



        % Gets entries where the defined parts of partial_key matches the key of the entry
        function [keys, entries] = find(obj,partial_key)
            keys = {};
            entries = {};
            for i = 1:length(obj.keys)
                if structIsSubset(partial_key,obj.keys{i})
                    i
                    obj.keys{i}
                    keys{end + 1}    = obj.keys{i};
                    entries{end + 1} = obj.getEntryByIndex(i);
                end
            end
        end
    end

    methods(Static)
        function merge(fn1, fn2, fnNew)
            sf1 = SolutionFile(fn1);
            sf2 = SolutionFile(fn2);

            sfNew = SolutionFile(fnNew);

            sfNew.keys = sf1.keys;
            sfNew.matfile.keys = sf1.keys;
            sfNew.matfile.entries = sf1.matfile.entries;

            for i = 1:length(sf2.keys)
                if sfNew.isKey(sf2.keys{i})
                    warning('Key ''%s'' exists in both files!',struct2string(sf2.keys{i}));
                end
                sfNew.store(sf2.keys{i},sf2.getEntryByIndex(i));
            end
        end


        function b = keyIsEqual(key1,key2)
            b = isequal(key1, key2);
        end

        function b = keyIsIn(key,keys)
            b = false;
            for i = 1:length(keys)
                b = isequal(key, keys{i});
                if b
                    return
                end
            end
        end
    end

end
