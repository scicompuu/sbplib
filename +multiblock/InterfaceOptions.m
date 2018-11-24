% An object of class InterfaceOptions can be passed as argument to multiblock.Diffop
% to specify details of the interface couplings.
%
% An InterfaceOptions object is essentially a cell array of options,
% equipped with methods that make it easier to change options.
classdef InterfaceOptions < handle
    properties
        % optsCell -- nBlocks x nBlocks cell array (same size as grid.connections)
        %             Must have the same sparsity pattern as grid.connections
        optsCell
    end

    methods

        % grid          --  mutliblock.grid
        % intialOpts    --  cell array of interface options
        function obj = InterfaceOptions(grid, initialOpts)

            default_arg('initialOpts', []);

            % If no initialOpts are given, create empty options cell.
            % The cell matrix is non-empty where connections is non-empty.
            if isempty(initialOpts)
                opts = grid.connections;
                for i = 1:numel(grid.connections)
                    if ~isempty(grid.connections{i})
                        opts{i} = cell(1, 2);
                        opts{i}{1} = struct;
                        opts{i}{2} = struct;
                    end
                end

            % If no grid is given, assume that initialOpts is correct and use it.
            elseif isempty(grid)
                opts = initialOpts;

            % Check that grid.connections and initialOpts match, and then use initialOpts.
            else
                assert(numel(grid.connections) == numel(initialOpts),...
                     'InterfaceOptions: grid.connections and initialOpts do not match');
                opts = initialOpts;
            end
            obj.optsCell = opts;
        end

        % Returns the cell matrix that contains the options
        function opts = getOptions(obj)
            opts = obj.optsCell;
        end


        % Sets the option optStr to val, for the coupling berween blocks i and j
        % If i and j are omitted, all couplings get optStr = val.
        %
        % optStr  -- string
        % val     -- anything
        % i,j     -- integers (or empty)
        function setOption(obj, optStr, val, i ,j)
            default_arg('i',[]);
            default_arg('j',[]);

            opts = obj.optsCell;

            if isempty(i) && ~isempty(j)
                error('If i is empty, j must also be empty.');

            elseif isempty(j) && ~isempty(i)
                error('If j is empty, i must also be empty.');

            % If i and j are empty, set the option for all interfaces
            elseif isempty(i) && isempty(j)
                for k = 1:numel(opts)
                    if ~isempty(opts{k})
                        opts{k}{1} = setfield(opts{k}{1}, optStr, val);
                        opts{k}{2} = setfield(opts{k}{2}, optStr, val);
                    end
                end

            % Both i and j are nonempty, set property only for that interface
            else
                if ~isempty(opts{i,j})
                    opts{i,j}{1} = setfield(opts{i,j}{1}, optStr, val);
                    opts{i,j}{2} = setfield(opts{i,j}{2}, optStr, val);
                elseif ~isempty(opts{j,i})
                    opts{j,i}{1} = setfield(opts{j,i}{1}, optStr, val);
                    opts{j,i}{2} = setfield(opts{j,i}{2}, optStr, val);
                else
                    error(sprintf('Blocks %d and %d do not seem to be coupled',i,j) );
                end
            end

            obj.optsCell = opts;
        end


        % Merges with another InterfaceOptions-object.
        % Errors if there are merge conflicts.
        % TODO: merge with preference?
        function merge(obj, obj2)
            localOpts = obj.getOptions();
            remoteOpts = obj2.getOptions();

            assert( numel(localOpts) == numel(remoteOpts), ...
                    'multiblock.InterfaceOptions: The two InterfaceOptions do not have the same dimension.');

            for i = 1:numel(localOpts)
                if ~isempty(remoteOpts{i})
                    if isempty(localOpts{i})
                        error('multiblock.InterfaceOptions: The two InterfaceOptions must have the same interface connections');
                    else
                        for j = 1:2
                            remoteStruct = remoteOpts{i}{j};
                            localStruct = localOpts{i}{j};
                            localFields = fieldnames(localStruct);

                            % Assert the we don't have any identical field names, which would lead to overwriting
                            for k = 1:numel(localFields)
                                if isfield(remoteStruct, localFields{k})
                                    error('multiblock.InterfaceOptions: Cannot perform union of InterfaceOptions with common options');
                                end
                            end

                            % Take fields from remote and deal to local
                            remoteFields = fieldnames(remoteStruct);
                            for k = 1:numel(remoteFields)
                                name = remoteFields{k};
                                val = getfield(remoteStruct, name);
                                localStruct = setfield(localStruct, name, val);
                            end

                            localOpts{i}{j} = localStruct;
                        end
                    end
                end
            end

            obj.optsCell = localOpts;
        end


    end
end
