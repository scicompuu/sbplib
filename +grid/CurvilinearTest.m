function tests = CurvilinearTest()
    tests = functiontests(localfunctions);
end

function testMappingInputGridFunction(testCase)
    in = {
        {{1:10}, @(x) exp(x)},
        {{1:10,1:6}, @(x,y) [exp(x+y); exp(x-y)]},
        {{1:10,1:5,1:7}, @(x,y,z)[exp(x+y+z); exp(x-y-z); 2+x+y-z]},
    };

    out = {
        [10, 1];
        [10*6, 2];
        [10*5*7, 3];
    };


    % How to test this? Just make sure it runs without errors.

    for i = 1:length(in)
        g = grid.Curvilinear(in{i}{2},in{i}{1}{:});
        testCase.verifyEqual(size(g.coords),out{i});
    end
end

function testMappingInputComponentMatrix(testCase)
    in = {
        {{1:3}, [1 2 3]'},
        {{1:2, 1:3}, [1 2 3 4 5 6; 7 8 9 10 11 12]'},
    };

    for i = 1:length(in)
        g = grid.Curvilinear(in{i}{2},in{i}{1}{:});
        testCase.verifyEqual(g.coords,in{i}{2});
    end
end

function testMappingInputCellOfMatrix(testCase)

    in = {
        {{1:3}, {[1 2 3]'}},
        {{1:2, 1:3}, {[1 2 3; 4 5 6], [7 8 9; 10 11 12]}},
    };

    out = {
        [1 2 3]',
        [1 2 3 4 5 6; 7 8 9 10 11 12]',
    };

    for i = 1:length(in)
        g = grid.Curvilinear(in{i}{2},in{i}{1}{:});
        testCase.verifyEqual(g.coords,out{i});
    end
end

function testMappingInputCellOfVectors(testCase)
    in = {
        {{1:3}, {[1 2 3]'}},
        {{1:2, 1:3}, {[1 2 3 4 5 6]', [7 8 9 10 11 12]'}},
    };

    out = {
        [1 2 3]',
        [1 2 3 4 5 6; 7 8 9 10 11 12]',
    };
end

function testMappingInputError(testCase)
    testCase.verifyFail();
end

function testScaling(testCase)
    in = {{1:2, 1:3}, {[1 2 3 4 5 6]', [7 8 9 10 11 12]'}};
    g = grid.Curvilinear(in{2},in{1}{:});

    testCase.verifyError(@()g.scaling(),'grid:Curvilinear:NoScalingSet');

    g.logicalGrid.h = [2 1];
    testCase.verifyEqual(g.scaling(),[2 1]);
end

function testGetBoundaryNames(testCase)
    in = {
        {{1:10}, @(x) exp(x)},
        {{1:10,1:6}, @(x,y) [exp(x+y); exp(x-y)]},
        {{1:10,1:5,1:7}, @(x,y,z)[exp(x+y+z); exp(x-y-z); 2+x+y-z]},
    };

    out = {
        {'l', 'r'},
        {'w', 'e', 's', 'n'},
        {'w', 'e', 's', 'n', 'd', 'u'},
    };

    for i = 1:length(in)
        g = grid.Curvilinear(in{i}{2},in{i}{1}{:});
        testCase.verifyEqual(g.getBoundaryNames(), out{i});
    end
end

function testGetBoundary(testCase)
    grids = {
        {{1:10}, @(x) exp(x)},
        {{1:10,1:6}, @(x,y) [exp(x+y); exp(x-y)]},
        {{1:10,1:5,1:7}, @(x,y,z)[exp(x+y+z); exp(x-y-z); 2+x+y-z]},
    };

    boundaries = {
        {'l', 'r'},
        {'w', 'e', 's', 'n'},
        {'w', 'e', 's', 'n', 'd', 'u'},
    };


    for ig = 1:length(grids)
        g = grid.Curvilinear(grids{ig}{2},grids{ig}{1}{:});

        logicalGrid = grid.Cartesian(grids{ig}{1}{:});

        for ib = 1:length(boundaries{ig})

            logicalBoundary = logicalGrid.getBoundary(boundaries{ig}{ib});

            x = num2cell(logicalBoundary',2);
            expectedBoundary = grids{ig}{2}(x{:})';
            testCase.verifyEqual(g.getBoundary(boundaries{ig}{ib}), expectedBoundary);
        end
    end
end

