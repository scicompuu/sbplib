function [d_min, i1_min, j1_min, i2_min, j2_min] = min_h(X,Y)
    ni = size(X,1);
    nj = size(X,2);
    d_min = norm([X(1,1);Y(1,1)] - [X(ni,nj);Y(ni,nj)]);

    i1_min = 0;
    j1_min = 0;
    i2_min = 0;
    j2_min = 0;

    D = {[-1,-1],[0,-1],[1,-1],[1,0],[1,1],[0,1],[-1,1],[-1,0]};
    % D = {[0,-1],[1,0],[0,1],[-1,0]};

    for i = 1:ni
        for j = 1:nj
            p1 = [X(i,j); Y(i,j)];
            for k = 1:length(D)
                i2 = i+D{k}(1);
                j2 = j+D{k}(2);
                if i2 >= 1 && i2 <= ni && j2 >= 1 && j2 <= nj
                    p2 = [X(i2,j2); Y(i2,j2)];
                    d = norm(p2-p1);
                    if d < d_min;
                        d_min = d;
                        i1_min = i;
                        j1_min = j;
                        i2_min = i2;
                        j2_min = j2;
                    end
                end
            end
        end
    end
end
