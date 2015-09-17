function [d_max, i1_max, j1_max, i2_max, j2_max] = max_h(X,Y)
    ni = size(X,1);
    nj = size(X,2);
    d_max = 0;

    i1_max = 0;
    j1_max = 0;
    i2_max = 0;
    j2_max = 0;

    D = {[0,-1],[1,0],[0,1],[-1,0]};

    for i = 1:ni
        for j = 1:nj
            p1 = [X(i,j); Y(i,j)];
            for k = 1:length(D)
                i2 = i+D{k}(1);
                j2 = j+D{k}(2);
                if i2 >= 1 && i2 <= ni && j2 >= 1 && j2 <= nj
                    p2 = [X(i2,j2); Y(i2,j2)];
                    d = norm(p2-p1);
                    if d > d_max;
                        d_max = d;
                        i1_max = i;
                        j1_max = j;
                        i2_max = i2;
                        j2_max = j2;
                    end
                end
            end
        end
    end
end
