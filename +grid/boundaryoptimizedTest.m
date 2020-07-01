function tests = boundaryoptimizedTest()
    tests = functiontests(localfunctions);
end

function testErrorInvalidParam(testCase)
    in  = {
        %Invalid order
        {[10 10],{0,1},{0,2},3},
        %Invalid grid size
        {5, {0,1}, 4},
        {[10 5],{0,1},{0,2},4},
        {[10 5],{0,1},{0,2},6,'M'},
        %Invalid limits
        {10,{1},4},
        {[10,10],{0,1},{1},4},
        {[10,10],{1},{1,0},4},
        {10,{1,0},4},
        {[10, 5],{1,0},{0,-1},4},
    };

    for i = 1:length(in)
        testCase.verifyError(@()grid.boundaryoptimized(in{i}{:}),'grid:boundaryoptimized:InvalidParameter',sprintf('in(%d) = %s',i,toString(in{i})));
    end
end

function testErrorInvalidOption(testCase)
    in  = {
        {[8 8],{0,1},{0,2},4,'acrurate'},
    };

    for i = 1:length(in)
        testCase.verifyError(@()grid.boundaryoptimized(in{i}{:}),'grid:boundaryoptimized:InvalidOption',sprintf('in(%d) = %s',i,toString(in{i})));
    end
end

function testErrorNonMatchingParam(testCase)
    in  = {
        {[],{1},4},
        {[],{0,1},{0,1},4},
        {[5,5],{0,1},{0,1},{0,1},4},
        {[5,5,4],{0,1},{0,1},4,'accurate'}
        {[5,5,4],{0,1},{0,1},{0,1},4,4},
    };

    for i = 1:length(in)
        testCase.verifyError(@()grid.boundaryoptimized(in{i}{:}),'grid:boundaryoptimized:NonMatchingParameters',sprintf('in(%d) = %s',i,toString(in{i})));
    end
end

% Tests that the expected grid points are obtained for a boundary optimized grid with a 4th order
% accurate stencil and 8th order minimal stencil.
% The boundary grid point distance weights are taken from the D1Nonequidistant operators and
% grid spacing is calculated according to Mattsson et al 2018. The test uses minimal number of grid
% points required by the operators.
function testCompiles(testCase)
    
    %% 1D 4th order accurate stencil
    % Boundary weights, number of non-equidistantly spaced points for 4th order accurate stencil
    bw = [0.0000000000000e+00 6.8764546205559e-01 1.8022115125776e+00];
    n = length(bw)-1;
    xi_n = bw(end);

    % Grid points in x-direction.
    Lx = 1;
    mx = 8;
    hx_4 = Lx/(2*xi_n+(mx-2*n-1)); 
    
    bp_l = hx_4*bw;
    bp_r = Lx-flip(hx_4*bw);
    interior = [hx_4*(xi_n+1) hx_4*(xi_n+2)];
    x_4 = [bp_l interior bp_r];

    % Boundary weights, number of non-equidistantly spaced points for 8th order minimal stencil    
    bw = [0.0000000000000e+00, 4.9439570885261e-01, 1.4051531374839e+00];
    n = length(bw)-1;
    xi_n = bw(end);

    %% 2D 8th order minimal stencil
    % Grid points in x-direction.
    hx_8 = Lx/(2*xi_n+(mx-2*n-1)); 
    
    bp_l = hx_8*bw;
    bp_r = Lx-flip(hx_8*bw);
    interior = [hx_8*(xi_n+1) hx_8*(xi_n+2)];
    x_8 = [bp_l interior bp_r];

    % Grid points in y-direction.
    Ly = 2;
    my = 9;
    hy = Ly/(2*xi_n+(my-2*n-1));
    
    bp_l = hy*bw;
    bp_r = Ly-flip(hy*bw);
    interior = [hy*(xi_n+1) hy*(xi_n+2) hy*(xi_n+3)];
    y = [bp_l interior bp_r];

    in  = {
        {mx, {0,Lx},4},
        {[mx, my],{0,Lx},{0,Ly},8,'M'},
    };
    
    out = {
        {[x_4'],hx_4}
        {[kr(x_8',ones(size(y'))),kr(ones(size(x_8')),y')],[hx_8, hy]}
    };

    for i = 1:length(in)
        g = grid.boundaryoptimized(in{i}{:});
        testCase.verifyEqual(g.points(),out{i}{1},'AbsTol', 1e-14, 'RelTol', 1e-14);
        testCase.verifyEqual(g.scaling(),out{i}{2},'AbsTol', 1e-14, 'RelTol', 1e-14);
    end
end