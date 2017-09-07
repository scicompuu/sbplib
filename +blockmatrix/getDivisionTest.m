function tests = getDivisionTest()
    tests = functiontests(localfunctions);
end

function testError(testCase)
    cases = {
        magic(3),
        {[2 2 2];{1,2}},
        {[2 2 2];[1 2]},
        {[2; 2; 2], [1; 2]},
    };

    for i =1:length(cases)
        testCase.verifyError(@()blockmatrix.getDivision(cases{i}), 'blockmatrix:getDivision:NotABlockmatrix')
    end
end

function testGetDivision(testCase)
    cases = {
        {
            {},
            {[],[]};
        },
        {
            {
            [2 2; 2 1], [1; 2];
            [2 2], [1]
            },
            {[2 1], [2 1]}
        },
        {
            {
            [2 2; 2 1], [];
            [2 2], [1]
            },
            {[2 1], [2 1]}
        },
        {
            {
            [2 2; 2 1], [];
            [2 2], []
            },
            {[2 1], [2 0]}
        },
        {
            {
            [2 2; 2 1], [1; 2];
            [], []
            },
            {[2 0], [2 1]}
        },
        {
            {
            [2 2; 2 1];
            [2 2]
            },
            {[2 1], 2}
        },
        {
            {zeros(3,0)},
            {3, 0},
        },
        {
            {zeros(3,0), zeros(3,0)},
            {3, [0, 0]},
        },
        {
            {zeros(3,0); zeros(2,0)},
            {[3 2],0},
        },
    };

    for i = 1:length(cases)
        in = cases{i}{1};
        out = blockmatrix.getDivision(in);
        expected = cases{i}{2};
        testCase.verifyEqual(out, expected);
    end
end


