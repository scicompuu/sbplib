function tests = isDivisionTest()
    tests = functiontests(localfunctions);
end

function testIsDivision(testCase)
    cases = {
        {{[2 2 2],[1 2]} ,true},
        {{[1 2],[1 0]} ,false},
        {{[0 2],[1 1]} ,false},
        {{[1 2],[]} ,false},
        {{[1 2],[1]} ,true},
        {{[1 2],[1], [1 2 3]} ,false},
        {{[1 2 3]} ,false},
        {[1 2] ,false},
    };

    for i = 1:length(cases)
        in = cases{i}{1};
        out = blockmatrix.isDivision(in);
        expected = cases{i}{2};
        testCase.verifyEqual(out, expected, sprintf('Should return %d for %s', expected, toString(in)));
    end
end