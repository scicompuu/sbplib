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
    cases = {
        {{},{}},
        {{1},{}},
        {{},{1}},
        {{1},{2}},
        {{1, 2},{3, 4}},
        {{1; 2},{3; 4}},
    };

    horzCat = {
        {},
        {1},
        {1},
        {1,2},
        {1, 2, 3, 4},
        {1, 3; 2, 4},
    };

    vertCat = {
        {},
        {1},
        {1},
        {1; 2},
        {1, 2; 3, 4},
        {1; 2; 3; 4},
    };

    for i = 1:length(cases)
        A = Cell(cases{i}{1});
        B = Cell(cases{i}{2});

        C_horz = [A, B];
        C_vert = [A; B];

        testCase.verifyEqual(C_horz.data, horzCat{i});
        testCase.verifyEqual(C_vert.data, vertCat{i});

    end
end
