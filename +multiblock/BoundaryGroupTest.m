function tests = BoundaryGroupTest()
    tests = functiontests(localfunctions);
end

function testCreation(testCase)
    in = {{3,'n'},{2,'hoho'},{1,'s'}};

    blockIDs = [3 2 1];
    names = {'n', 'hoho', 's'};

    bg = multiblock.BoundaryGroup(in{:});
    testCase.verifyEqual(bg.blockIDs, blockIDs);
    testCase.verifyEqual(bg.names, names);
end

function testInputError(testCase)
    in = {
        {'n', 's'},
        {{3,'n'},{2,2,'hoho'},{1,'s'}},
    };

    out = {
        'multiblock:BoundaryGroup:BoundaryGroup:InvalidInput',
        'multiblock:BoundaryGroup:BoundaryGroup:InvalidInput',
    };

    for i = 1:length(in)
        testCase.verifyError(@()multiblock.BoundaryGroup(in{i}{:}), out{i});
    end
end