% Stitch schemes together given connection matrix and BC vector.
%     schmHand  - function_handle to a Scheme constructor
%     order     - order of accuracy
%     schmParam - cell array of extra parameters sent to each Scheme stored as cell arrays
%     blocks    - block definitions, On whatever form the scheme expects.
%     ms        - grid points in each direction for each block. Ex {[10,10], [10, 20]}
%     conn      - connection matrix
%     bound     - boundary condition vector, array of structs with fields w,e,s,n
%                 each field with a parameter array that is sent to schm.boundary_condition
%
% Output parameters are cell arrays and cell matrices.
%
% Ex: [schms, D, H] = stitchSchemes(schmHand, order, schmParam, blocks, ms, conn, bound)
function [schms, D, H] = stitchSchemes(schmHand, order, schmParam, blocks, ms, conn, bound)
    default_arg('schmParam',[]);

    n_blocks = numel(blocks);

    % Creating Schemes
    for i = 1:n_blocks
        if isempty(schmParam);
            schms{i} = schmHand(ms{i},blocks{i},order,[]);
        elseif ~iscell(schmParam)
            param = schmParam(i);
            schms{i} = schmHand(ms{i},blocks{i},order,param);
        else
            param = schmParam{i};
            if iscell(param)
                schms{i} = schmHand(ms{i},blocks{i},order,param{:});
            else
                schms{i} = schmHand(ms{i},blocks{i},order,param);
            end
        end

        % class(schmParam)
        % class(ms)
        % class(blocks)
        % class(schmParam{i})
        % class(ms)


    end


    % Total norm
    H = cell(n_blocks,n_blocks);
    for i = 1:n_blocks
        H{i,i} = schms{i}.H;
    end

    %% Total system matrix

    % Differentiation terms
    D = cell(n_blocks,n_blocks);
    for i = 1:n_blocks
        D{i,i} = schms{i}.D;
    end

    % Boundary penalty terms
    for i = 1:n_blocks
        if ~isstruct(bound{i})
            continue
        end

        fn = fieldnames(bound{i});
        for j = 1:length(fn);
            bc = bound{i}.(fn{j});
            if isempty(bc)
                continue
            end

            t = schms{i}.boundary_condition(fn{j},bc{:});
            D{i,i} = D{i,i}+t;
        end
    end

    % Interface penalty terms
    for i = 1:n_blocks
        for j = 1:n_blocks
            intf = conn{i,j};
            if isempty(intf)
                continue
            end

            [uu,uv,vv,vu] = schms{i}.interface_coupling(schms{i},intf{1},schms{j},intf{2});
            D{i,i} = D{i,i} + uu;
            D{i,j} = uv;
            D{j,j} = D{j,j} + vv;
            D{j,i} = vu;
        end
    end
end