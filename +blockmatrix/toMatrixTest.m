function tests = toMatrixTest()
    tests = functiontests(localfunctions);
end

function testError(testCase)
    testCase.verifyError(@()blockmatrix.toMatrix([]), 'blockmatrix:toMatrix:NotABlockmatrix')
end

function testToMatrix(testCase)
    cases = {
        {
            {},
            [],
        },
        {
            {1, 2; 3, 4},
            [1,2; 3,4],
        }
        {
            {
            [2 2; 2 1], [1; 2];
            [2 2], [1]
            },
            [2 2 1;
             2 1 2;
             2 2 1],
        },
        {
            {
            [2 2; 2 1], [];
            [2 2], [1]
            },
            [2 2 0;
             2 1 0;
             2 2 1],
        },
        {
            {
            [2 2; 2 1], [];
            [2 2], []
            },
            [2 2;
             2 1;
             2 2],
        },
        {
            {
            [2 2; 2 1], [1; 2];
            [], []
            },
            [2 2 1;
             2 1 2],
        },
        {
            {zeros(0,0)},
            [],
        },
        {
            {zeros(3,0), zeros(3,0)},
            zeros(3,0),
        },
        {
            {zeros(3,0); zeros(2,0)},
            zeros(5,0),
        },
        {
            {zeros(0,3), zeros(0,2)},
            zeros(0,5),
        },
    };

    for i = 1:length(cases)
        in = cases{i}{1};
        out = full(blockmatrix.toMatrix(in));
        expected = cases{i}{2};
        testCase.verifyEqual(out, expected);
    end
end
