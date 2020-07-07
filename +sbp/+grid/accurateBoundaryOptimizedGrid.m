% Computes the grid points x and grid spacing h used by the boundary optimized SBP operators 
% with improved boundary accuracy, presented in 
% 'Boundary optimized diagonal-norm SBP operators - Mattsson, Almquist, van der Weide 2018'.
%
% lim - cell array with domain limits
% N - Number of grid points
% order - order of accuracy of sbp operator.
function [x,h] = accurateBoundaryOptimizedGrid(lim,N,order)
    assert(iscell(lim) && numel(lim) == 2,'The limit should be cell array with 2 elements.');
    L = lim{2} - lim{1};
    assert(L>0,'Limits must be given in increasing order.');
    %%%% Non-equidistant grid points %%%%%
    xb = boundaryPoints(order);
    m = length(xb)-1; % Number of non-equidistant points
    assert(N-2*(m+1)>=0,'Not enough grid points to contain the boundary region. Requires at least %d points.',2*(m+1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% Compute h %%%%%%%%%%
    h = L/(2*xb(end) + N-1-2*m);
    %%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% Define grid %%%%%%%%
    x = h*[xb; linspace(xb(end)+1,L/h-xb(end)-1,N-2*(m+1))'; L/h-flip(xb) ];
    x = x + lim{1};
    %%%%%%%%%%%%%%%%%%%%%%%%%
end
function xb = boundaryPoints(order)
    switch order
        case 4
            x0 =  0.0000000000000e+00;
            x1 =  6.8764546205559e-01;
            x2 =  1.8022115125776e+00;
            xb = [x0 x1 x2]';
        case 6
            x0 =  0.0000000000000e+00;
            x1 =  4.4090263368623e-01;
            x2 =  1.2855984345073e+00;
            x3 =  2.2638953951239e+00;
            xb = [x0 x1 x2 x3]';
        case 8
            x0 =  0.0000000000000e+00;
            x1 =  3.8118550247622e-01;
            x2 =  1.1899550868338e+00;
            x3 =  2.2476300175641e+00;
            x4 =  3.3192851303204e+00;
            xb = [x0 x1 x2 x3 x4]';
        case 10
            x0 =  0.0000000000000e+00;
            x1 =  3.5902433622052e-01;
            x2 =  1.1436659188355e+00;
            x3 =  2.2144895894456e+00;
            x4 =  3.3682742337736e+00;
            x5 =  4.4309689056870e+00;
            xb = [x0 x1 x2 x3 x4 x5]';
        case 12
            x0 =  0.0000000000000e+00;
            x1 =  3.6098032343909e-01;
            x2 =  1.1634317168086e+00;
            x3 =  2.2975905356987e+00;
            x4 =  3.6057529790929e+00;
            x5 =  4.8918275675510e+00;
            x6 =  6.0000000000000e+00;
            xb = [x0 x1 x2 x3 x4 x5 x6]';
        otherwise
            error('Invalid operator order %d.',order);
    end
end