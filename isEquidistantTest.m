function tests = isEquidistantTest()
    tests = functiontests(localfunctions);
end

function testTooShortInput(testCase)
    testCase.verifyError(@()isEquidistant([]), 'sbplib:isEquidistant:inputTooShort')
end

function testCorrectOutput(testCase)
    cases = {
        % {input, expected},
        {[0,0,0,0,0], true},
        {[1,1,1,1,1], true},
        {[1,2,3,4,5], true},
        {[1,3,4,5], false},
        {[1,2,3,5], false},
        {[1,2,4,5], false},
        {linspace(0,pi, 3), true},
        {linspace(0,1, 4), true},
        {linspace(0,1, 4123), true},
        {linspace(0,sin(1), 123), true},
    };

    for i = 1:length(cases)
        input = cases{i}{1};
        expected = cases{i}{2};
        result = isEquidistant(input);

        testCase.verifyEqual(result,expected);
    end
end
