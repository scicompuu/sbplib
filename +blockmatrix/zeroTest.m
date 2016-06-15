function tests = zeroTest()
    tests = functiontests(localfunctions);
end

function testZero(testCase)
    cases = {
        {
            {[],[]},
            {},
        },
        {
            {0,0},
            {[]};
        },
        {
            {1,1},
            {0};
        },
        {
            {2,1},
            {[0; 0]};
        },
        {
            {1,2},
            {[0 0]};
        },
        {
            {[1 2],2},
            {[0 0];[0 0; 0 0]};
        },
        {
            {[1 2],[2 1]},
            {[0 0],[0];[0 0; 0 0],[0; 0]};
        },
    };

    for i = 1:length(cases)
        out = convertToFull(blockmatrix.zero(cases{i}{1}));
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