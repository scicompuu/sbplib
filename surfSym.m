function h = surfSym(X,Y,Z)
    if isvector(X) && isvector(Y)
        [X,Y] = meshgrid(X,Y);
    end

    V_nodes = [X(:), Y(:), Z(:)];

    X_centers = neighbourMean(X);
    Y_centers = neighbourMean(Y);
    Z_centers = neighbourMean(Z);
    V_centers = [X_centers(:), Y_centers(:), Z_centers(:)];


    N = prod(size(X));
    nodeIndecies = reshape(1:N, size(X));
    centerIndecies = reshape(N+(1:prod(size(X)-[1,1])), size(X) - [1,1]);

    % figure()
    % h = line(V_nodes(:,1),V_nodes(:,2),V_nodes(:,3));
    % h.LineStyle = 'none';
    % h.Marker = '.';
    % h = line(V_centers(:,1),V_centers(:,2),V_centers(:,3));
    % h.LineStyle = 'none';
    % h.Marker = '.';
    % h.Color = Color.red;
    % axis equal


    I_0 = nodeIndecies(1:end-1, 1:end-1);
    I_1 = nodeIndecies(2:end, 1:end-1);
    I_2 = nodeIndecies(2:end, 2:end);
    I_3 = nodeIndecies(1:end-1, 2:end);
    I_C = centerIndecies;

    S.Vertices = [
        V_nodes;
        V_centers;
    ];

    S.Faces = [
        I_0(:), I_1(:), I_C(:);
        I_1(:), I_2(:), I_C(:);
        I_2(:), I_3(:), I_C(:);
        I_3(:), I_0(:), I_C(:);
    ];

    % figure()
    h = patch(S, 'FaceVertexCData',  S.Vertices(:,3),'FaceColor','flat');
end

% Calculate the mean of four neighbours around a patch
function M = neighbourMean(A)
    M = (A(1:end-1, 1:end-1) + A(2:end, 1:end-1) + A(1:end-1, 2:end) + A(2:end, 2:end))/4;
end
