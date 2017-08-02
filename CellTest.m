function tests = CellTest()
    tests = functiontests(localfunctions);
end

function testSubAssignment(testCase)
    testCase.verifyFail();
end

function testIndexreference(testCase)
    testCase.verifyFail();
end

function testConcat(testCase)
    testCase.verifyFail();
end
