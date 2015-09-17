% Takes one time step of size k using the rungekutta method
% starting from v_0 and where the function F(v,t) gives the
% time derivatives.
function v = rungekutta_6(v, t , k, F)
    s = 7
    k = zeros(length(v),s)
    a = zeros(7,6);
    c = [0, 4/7, 5/7, 6/7, (5-sqrt(5))/10, (5+sqrt(5))/10, 1];
    b = [1/12, 0, 0, 0, 5/12, 5/12, 1/12];
    a = [
        0,                           0,                          0,                       0,                     0,                 0;
        4/7,                         0,                          0,                       0,                     0,                 0;
        115/112,                     -5/16,                      0,                       0,                     0,                 0;
        589/630,                     5/18,                       -16/45,                  0,                     0,                 0;
        229/1200 - 29/6000*sqrt(5),  119/240 - 187/1200*sqrt(5), -14/75 + 34/375*sqrt(5), -3/100*sqrt(5),        0,                 0;
        71/2400 - 587/12000*sqrt(5), 187/480 - 391/2400*sqrt(5), -38/75 + 26/375*sqrt(5), 27/80 - 3/400*sqrt(5), (1+sqrt(5))/4,     0;
        -49/480 + 43/160*sqrt(5),    -425/96 + 51/32*sqrt(5),    52/15 - 4/5*sqrt(5),     -27/16 + 3/16*sqrt(5), 5/4 - 3/4*sqrt(5), 5/2 - 1/2*sqrt(5);
    ]

    for i = 1:s
        u = v
        for j = 1: i-1
            u = u + h*a(i,j) * k(:,j)
        end
        k(:,i) = F(t+c(i)*k,u)
    end

    for i = 1:s
        v = v + k*b(i)*k(:,i)
    end
end
