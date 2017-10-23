function tests = TiTest()
    tests = functiontests(localfunctions);
end

function testScalarInput(testCase)
    ti = getMinimumTi();

    cases = {
        % {u, v, out},
        {0, 0, [1; 2]},
        {0, 1, [1; 4]},
        {1, 0, [3; 2]},
        {1, 1, [3; 4]},
        {0.5, 0.5, [2; 3]},
    };

    for i = 1:length(cases)
        u = cases{i}{1};
        v = cases{i}{2};
        expected = cases{i}{3};

        testCase.verifyEqual(ti.S(u,v), expected, sprintf('Case: %d',i));
    end
end

function testRowVectorInput(testCase)
    ti = getMinimumTi();

    u = [0, 0.5, 1];
    v = [0, 0, 0.5];
    expected = [
        1, 2, 3;
        2, 2, 3;
    ];

    testCase.verifyEqual(ti.S(u,v), expected);
end

function testColumnvectorInput(testCase)
   ti = getMinimumTi();

    u = [0; 0.5; 1];
    v = [0; 0; 0.5];
    expected = [1; 2; 3; 2; 2; 3];

    testCase.verifyEqual(ti.S(u,v), expected);
end


function ti = getMinimumTi()
    ti = parametrization.Ti.rectangle([1; 2], [3; 4]);
end