function tests = equdistantCurvilinearTest()
    tests = functiontests(localfunctions);
end

function testNoTests(testCase)
    testCase.assumeFail();
end