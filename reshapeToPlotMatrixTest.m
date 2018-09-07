function tests = reshapeToPlotMatrixTest()
    tests = functiontests(localfunctions);
end

function test1D(testCase)
    inGf = [1 2 3 4 5]';
    inM = 5;
    out = [1 2 3 4 5]';
    testCase.verifyEqual(reshapeToPlotMatrix(inGf, inM),out);
end

function test2D(testCase)
    x = 1:2;
    y = 1:3;

    f = @(x,y) x + y*10;

    xx = [1; 1; 1; 2; 2; 2];
    yy = [1; 2; 3; 1; 2; 3];
    inGf = f(xx,yy);

    [X,Y] = meshgrid(x,y);
    out = f(X,Y);

    inM = [2, 3];

    testCase.verifyEqual(reshapeToPlotMatrix(inGf, inM),out);
end

function test3D(testCase)
    x = 1:2;
    y = 1:3;
    z = 1:4;

    f = @(x,y,z) x + y*10 + z*100;

    xx = [repmat(1, [12, 1]); repmat(2, [12, 1])];
    yy = repmat([1; 1; 1; 1; 2; 2; 2; 2; 3; 3; 3; 3], [2, 1]);
    zz = repmat([1; 2; 3; 4], [6, 1]);
    inGf = f(xx,yy,zz);

    [X,Y,Z] = meshgrid(x,y,z);
    out = f(X,Y,Z);

    inM = [2, 3, 4];

    testCase.verifyEqual(reshapeToPlotMatrix(inGf, inM),out);
end