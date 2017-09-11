function tests = isDivisionTest()
    tests = functiontests(localfunctions);
end

function testIsDivision(testCase)
    cases = {
        {[1 2] ,false},     % Must be a cell array
        {{[1 2 3]} ,false}, % Must have two vectors
        {{[],[]}, true}     % No blocks is a valid blockmatrix
        {{[1 2],[]} ,true},
        {{[],[1 2]} ,true},
        {{[2 2 2],[1 2]} ,true},
        {{[1 2],[1 0]} ,true},
        {{[0 2],[1 1]} ,true},
        {{[1 2],[1]} ,true},
        {{[1 2],[1], [1 2 3]} ,false},
    };

    for i = 1:length(cases)
        in = cases{i}{1};
        out = blockmatrix.isDivision(in);
        expected = cases{i}{2};
        testCase.verifyEqual(out, expected, sprintf('Should return %d for %s', expected, toString(in)));
    end
end