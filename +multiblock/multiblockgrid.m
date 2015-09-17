% Creates a multi block square grid with defined boundary conditions.
%   x,y defines the grid lines. Rember to think of the indexing as a matrix. Order matters!
%   bc is a struct defining the boundary conditions on each side of the block.
%       bc.w = {'dn',[function or value]}
function [block,conn,bound,ms] = multiblockgrid(x,y,mx,my,bc)
    n = length(y)-1; % number of blocks in the y direction.
    m = length(x)-1; % number of blocks in the x direction.
    N = n*m; % number of blocks

    if ~issorted(x)
        error('The elements of x seem to be in the wrong order');
    end
    if ~issorted(flip(y))
        error('The elements of y seem to be in the wrong order');
    end
    % y = sort(y,'descend');

    % Dimensions of blocks and number of points
    block = cell(n,m);
    for i = 1:n
        for j = 1:m
            block{i,j} = {
                {x(j),x(j+1)}, {y(i+1),y(i)};
            };

            ms{i,j} = [mx(i),my(j)];
        end
    end

    % Interface couplings
    conn = cell(N,N);
    for i = 1:n
        for j = 1:m
            I = flat_index(n,i,j);
            if i < n
                J = flat_index(n,i+1,j);
                conn{I,J} = {'s','n'};
            end

            if j < m
                J = flat_index(n,i,j+1);
                conn{I,J} = {'e','w'};
            end
        end
    end


    % Boundary conditions
    bound = cell(n,m);
    for i = 1:n
        if isfield(bc,'w')
            bound{i,1}.w = bc.w;
        end

        if isfield(bc,'e')
            bound{i,n}.e = bc.e;
        end
    end

    for j = 1:m
        if isfield(bc,'n')
            bound{1,j}.n = bc.n;
        end

        if isfield(bc,'s')
            bound{m,j}.s = bc.s;
        end
    end
end

