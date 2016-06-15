function tests = fromMatrixTest()
    tests = functiontests(localfunctions);
end

function testErrorNonMatchingDim(testCase)
    in  = {
        {magic(5), {[1 2 3], [4]}},
        {magic(5), {[1 1 1 1 1 1], [5]}},
        {magic(5), {[5], [1 1 1 1 1 1]}},
        {ones(4,2),{[2 3],[2]}},
        {ones(4,2),{[2 2],[3]}},
    };

    for i = 1:length(in)
        testCase.verifyError(@()blockmatrix.fromMatrix(in{i}{:}),'blockmatrix:fromMatrix:NonMatchingDim');
    end
end

function testFromMatrix(testCase)
    cases = {
        {
            {[],{[],[]}},
            {}
        },
        {
            {
                magic(3),
                {[3],[3]}
            },
            {magic(3)}
        },
        {
            {
                magic(3),
                {[1 1 1],[1 1 1]}
            },
            mat2cell(magic(3),[1 1 1],[1 1 1])
        },
        {
            {
                [17 24 1 8 15; 23 5 7 14 16; 4 6 13 20 22; 10 12 19 21 3; 11 18 25 2 9],
                {[1 4],[2 3]}
            },
            {
                [17 24], [1 8 15];
                [23 5; 4 6; 10 12; 11 18], [7 14 16; 13 20 22; 19 21 3; 25 2 9];
            };
        },
    };
    for i = 1:length(cases)
        out = convertToFull(blockmatrix.fromMatrix(cases{i}{1}{:}));
        expected = cases{i}{2};
        testCase.verifyEqual(out,expected);
    end
end

function C = convertToFull(C)
    [N,M] = size(C);
    for i = 1:N
        for j = 1:M
            C{i,j} = full(C{i,j});
        end
    end
end