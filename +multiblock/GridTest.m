function tests = GridTest()
    tests = functiontests(localfunctions);
end

function testCreation(testCase)
    g = multiblock.Grid({},{});
end

function testMissing(testCase)
    testCase.verifyFail();
end

function testGetBoundaryNames(testCase)
    testCase.verifyFail();
end

function testGetBoundary(testCase)
    testCase.verifyFail();
end