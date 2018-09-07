function tests = evalOnScalarTest()
    tests = functiontests(localfunctions);
end

function testInputConstant(testCase)
    in  = {
        0,
        47,
        1,
        [1; 2],
    };

    out = {
        [0; 0; 0],
        [47; 47; 47],
        [1; 1; 1],
        [1; 2; 1; 2; 1; 2],
    };

    g = getTestGrid('1d');

    for i = 1:length(in)
        gf = grid.evalOnScalar(g,in{i});
        testCase.verifyEqual(gf, out{i});
    end
end

function testInputScalarFunction1d(testCase)
    in  = {
        @(x)1+x*0,
        @(x)x,
        @(x)x.*x,
    };

    out = {
        [1; 1; 1],
        [0; 1; 2],
        [0; 1; 4],
    };

    g = getTestGrid('1d');

    for i = 1:length(in)
        gf = grid.evalOnScalar(g,in{i});
        testCase.verifyEqual(gf, out{i});
    end
end

function testInputScalarFunction2d(testCase)
    in  = {
        @(x,y)1+x*0,
        @(x,y)x-y,
        @(x,y)x./(1+y),
    };

    out = {
        [1; 1; 1; 1; 1; 1; 1; 1; 1],
        [0; -1; -2; 1; 0; -1; 2; 1; 0],
        [0; 0; 0; 1; 1/2; 1/3; 2; 1; 2/3],
    };

    g = getTestGrid('2d');

    for i = 1:length(in)
        gf = grid.evalOnScalar(g, in{i});
        testCase.verifyEqual(gf, out{i});
    end
end


function testInputVectorFunction(testCase)
    g = getTestGrid('1d');
    in = @(x)[x; -2*x];
    out = [0; 0; 1; -2; 2; -4];

    gf = grid.evalOnScalar(g,in);
    testCase.verifyEqual(gf, out);

    g = getTestGrid('2d');
    in = @(x,y)[x.^2; -2*y];
    out = [
        0;  0;
        0; -2;
        0; -4;
        1;  0;
        1; -2;
        1; -4;
        4;  0;
        4; -2;
        4; -4;
    ];

    gf = grid.evalOnScalar(g,in);
    testCase.verifyEqual(gf, out);
end


function testInputErrorVectorValued(testCase)
     in  = {
        [1,2,3],
        @(x,y)[x,-y];
    };

    g = getTestGrid('2d');

    for i = 1:length(in)
        testCase.verifyError(@()grid.evalOnScalar(g, in{i}),'grid:evalOnScalar:VectorValuedWrongDim',sprintf('in(%d) = %s',i,toString(in{i})));
    end
end

function g = getTestGrid(d)
    switch d
        case '1d'
            g = grid.equidistant(3,{0,2});
        case '2d'
            g = grid.equidistant([3,3],{0,2},{0,2});
    end
end
