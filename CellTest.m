function tests = CellTest()
    tests = functiontests(localfunctions);
end

function testSize(testCase)
    cases = {
        {{}, [0, 0]},
        {{1}, [1, 1]},
        {{1, 2}, [1, 2]},
        {{1; 2}, [2, 1]},
        {{1, 2; 3, 4}, [2,2]},
    };

    for i = 1:length(cases)
        A = Cell(cases{i}{1});
        expected = cases{i}{2};

        testCase.verifyEqual(size(A),expected);
    end
end

function testLength(testCase)
    cases = {
        {{}, 0},
        {{1}, 1},
        {{1, 2}, 2},
        {{1; 2}, 2},
        {{1, 2; 3, 4}, 2},
    };

    for i = 1:length(cases)
        A = Cell(cases{i}{1});
        expected = cases{i}{2};

        testCase.verifyEqual(length(A),expected);
    end
end

function testIsEmpty(testCase)
    cases = {
        {cell(0,0), true},
        {cell(1,0), true},
        {cell(0,1), true},
        {cell(1,1), false},
    };

    for i = 1:length(cases)
        A = Cell(cases{i}{1});
        expected = cases{i}{2};
        testCase.verifyEqual(isempty(A),expected);
    end
end

function testTranspose(testCase)
    testCase.verifyEqual(Cell({1i, 2}).', Cell({1i; 2}));
    testCase.verifyEqual(Cell({1i; 2}).', Cell({1i, 2}));
end

function testCtranspose(testCase)
    testCase.verifyEqual(Cell({1i, 2})', Cell({1i; 2}));
    testCase.verifyEqual(Cell({1i; 2})', Cell({1i, 2}));
end

function testRoundIndexWithProperty(testCase)
    A = Cell({3,2,1});

    result = A([1,3]).data;
    testCase.verifyEqual(result, {3, 1});
end

function testSubAssignmentRound(testCase)
    cases = {
        % {
        %     lArray,
        %     index,
        %     rhs,
        %     expectedResult
        % },
        {
            {},
            1,
            {'a'},
            {'a'},
        },
        {
            {1},
            1,
            {'a'},
            {'a'},
        },
        {
            {1,2,3},
            2,
            {'a'},
            {1,'a',3},
        },
        {
            {1,2,3},
            2,
            [],
            {1,3},
        },
    };

    for i = 1:length(cases)
        lArray         = Cell(cases{i}{1});
        index          = cases{i}{2};
        rhs            = cases{i}{3};
        expectedResult = cases{i}{4};

        lArray(index) = rhs;

        testCase.verifyEqual(lArray.data, expectedResult)
    end
end

function testSubAssignmentCurly(testCase)
    cases = {
        % {
        %     lArray,
        %     index,
        %     rhs,
        %     expectedResult
        % },
        {
            {},
            1,
            'a',
            {'a'},
        },
        {
            {1},
            1,
            'a',
            {'a'},
        },
        {
            {1,2,3},
            2,
            'a',
            {1,'a',3},
        },
    };

    for i = 1:length(cases)
        lArray         = Cell(cases{i}{1});
        index          = cases{i}{2};
        rhs            = cases{i}{3};
        expectedResult = cases{i}{4};

        lArray{index} = rhs;

        testCase.verifyEqual(lArray.data, expectedResult)
    end
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

    testCase.verifyEqual(C(:), D);
    testCase.verifyEqual(D(:), D);
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
