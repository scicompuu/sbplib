function tests = funcToComponentsTest()
    tests = functiontests(localfunctions);
end


function testScalarGf(testCase)
    g = getTestGrid();
    gf_in = [1; 2; 3];

    testCase.verifyEqual(grid.funcToComponents(g, gf_in), gf_in);
end

function testVectorGf(testCase)
    g = getTestGrid();
    gf_in = [1; 2; 3; 4; 5; 6];
    out = [1 2; 3 4; 5 6];

    testCase.verifyEqual(grid.funcToComponents(g, gf_in), out);
end

function g = getTestGrid()
    g = grid.equidistant(3,{0,2});
end