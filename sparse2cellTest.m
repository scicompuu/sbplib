function tests = sparse2cellTest()
    tests = functiontests(localfunctions);
end

function testErrorNonMatchingDim(testCase)
    in  = {
        {magic(5), [1 2 3], [4]},
        {magic(5), [1 1 1 1 1 1], [5]},
        {magic(5), [5], [1 1 1 1 1 1]},
        {ones(4,2),[2 3],[2]},
        {ones(4,2),[2 2],[3]},
    };

    for i = 1:length(in)
        testCase.verifyError(@()sparse2cell(in{i}{:}),'sparse2cell:NonMatchingDim');
    end
end

function testOutput(testCase)
    in = {};
    out = {};
    in{1}{1} =[17 24 1 8 15; 23 5 7 14 16; 4 6 13 20 22; 10 12 19 21 3; 11 18 25 2 9];
    in{1}{2} = [1 4];
    in{1}{3} = [2 3];

    out{1} = {
        [17 24], [1 8 15];
        [23 5; 4 6; 10 12; 11 18], [7 14 16; 13 20 22; 19 21 3; 25 2 9];
    }

    in{1}{1} = [17 24 1 8 15; 23 5 0 0 0; 4 6 0 0 0; 10 12 0 0 0; 11 18 0 0 0];
    in{1}{2} = [1 4];
    in{1}{3} = [2 3];

    out{1} = {
        [17 24], [1 8 15];
        [23 5; 4 6; 10 12; 11 18], [];
    }

    for i = 1:length(in)
        testCase.verifyEqual(sparse2cell(in{i}{:}), out{i});
    end
end