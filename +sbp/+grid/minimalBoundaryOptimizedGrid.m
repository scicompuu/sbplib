% Computes the grid points x and grid spacing h used by the boundary optimized SBP operators 
% with minimal number of non-equidistant boundary points, presented in 
% 'Boundary optimized diagonal-norm SBP operators - Mattsson, Almquist, van der Weide 2018'.
%
% lim - cell array with domain limits
% N - Number of grid points
% order - order of accuracy of sbp operator.
function [x,h] = minimalBoundaryOptimizedGrid(lim,N,order)
	assert(iscell(lim) && numel(lim) == 2,'The limit should be a cell array with 2 elements.');
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
			x1 =  7.7122987842562e-01;
			xb = [x0 x1]';
		case 6
			x0 =  0.0000000000000e+00;
			x1 =  4.0842950991998e-01;
			x2 =  1.1968523189207e+00;
			xb = [x0 x1 x2]';
		case 8
			x0 =  0.0000000000000e+00;
			x1 =  4.9439570885261e-01;
			x2 =  1.4051531374839e+00;
			xb = [x0 x1 x2]';
		case 10
			x0 =  0.0000000000000e+00;
			x1 =  5.8556160757529e-01;
			x2 =  1.7473267488572e+00;
			x3 =  3.0000000000000e+00;
			xb = [x0 x1 x2 x3]';
		case 12
			x0 =  0.0000000000000e+00;
			x1 =  4.6552112904489e-01;
			x2 =  1.4647984306493e+00;
			x3 =  2.7620429464763e+00;
			x4 =  4.0000000000000e+00;
			xb = [x0 x1 x2 x3 x4]';
		otherwise
			error('Invalid operator order %d.',order);
	end
end