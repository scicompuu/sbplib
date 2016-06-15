function tests = isBlockmatrixTest()
    tests = functiontests(localfunctions);
end

function testIsBlockmatrix(testCase)
    cases = {
        {
            magic(3),
            false % Must be a cell array
        }
        {
            {[2 2 2];{1,2}},
            false % All elements of the cell matrix must be regular matrices
        },
        {
            {[2 2 2];[1 2]},
            false % Row dimensions must match
        },
        {
            {[2; 2; 2], [1; 2]},
            false % Column dimensions must match
        },
        {
            {
            [2 2; 2 1], [1; 2];
            [2 2], [1]
            },
            true % A simple valid one
        },
        {
            {
            [2 2; 2 1], [];
            [2 2], [1]
            },
            true % Empty blocks assumed to be zero and match dimensions
        },
        {
            {
            [2 2; 2 1], [];
            [2 2], []
            },
            true % Empty blocks allowed.
        },
        {
            {
            [2 2; 2 1], [1; 2];
            [], []
            },
            true % Empty blocks allowed.
        },


    };

    for i = 1:length(cases)
        in = cases{i}{1};
        out = blockmatrix.isBlockmatrix(in);
        expected = cases{i}{2};
        testCase.verifyEqual(out, expected, sprintf('Should return %d for %s', expected, toString(in)));
    end
end