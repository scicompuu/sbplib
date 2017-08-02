function tests = CellTest()
    tests = functiontests(localfunctions);
end

function testSize(testCase)
    testCase.verifyFail();
end

function testLength(testCase)
    testCase.verifyFail();
end

function testTranspose(testCase)
    testCase.verifyFail();
end

function testRoundIndexWithProperty(testCase)
    A = Cell({3,2,1});

    testCase.verifyEqual(A([1,3]).data, {3, 1});
end

function testSubAssignment(testCase)
    testCase.verifyFail();
end

function testIndexreferenceRound(testCase)
    cases = {
        % {
        %     array,
        %     index,
        %     roundResult
        % },
        {
            {1,2,3},
            1,
            {1},
        },
        {
            {1,3,2},
            2,
            {3},
        },
        {
            {1,3,2},
            [1 3],
            {1, 2},
        },
    };


    for i = 1:length(cases)
        array = Cell(cases{i}{1});
        index = cases{i}{2};
        expected = cases{i}{3};

        result = array(index);

        testCase.verifyTrue(isa(result, 'Cell'));
        testCase.verifyEqual(result.data, expected);
    end
end

function testEndIndexing(testCase)
    C = Cell({1,2,3});

    testCase.verifyEqual(C(end), Cell({3}));
    testCase.verifyEqual(C{end}, 3);
end

function testColonIndexing(testCase)
    C = Cell({1, 2, 3});
    D = Cell({1; 2; 3});

    testCase.verifyEqual(C(:), Cell({3}));


    testCase.verifyEqual(C(:), Cell({3}));
    testCase.verifyEqual(C{end}, 3);
end

function testIndexreferenceCurly(testCase)
    cases = {
        % {
        %     array,
        %     index,
        %     curlyResult
        % },
        {
            {1,2,3},
            1,
            1
        },
        {
            {1,3,2},
            2,
            3
        },
    };


    for i = 1:length(cases)
        array = Cell(cases{i}{1});
        index = cases{i}{2};
        expected = cases{i}{3};

        result = array{index};

        testCase.verifyEqual(result, expected);
    end
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
