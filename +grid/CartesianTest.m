function tests = CartesianTest()
    tests = functiontests(localfunctions);
end


function testWarningEmptyGrid(testCase)
    in  = {
        {[]},
        {[],[1]},
        {[1],[2], []},
    };

    for i = 1:length(in)
        testCase.verifyError(@()grid.Cartesian(in{i}{:}),'grid:Cartesian:EmptyGrid');
    end
end

function testN(testCase)
    in  = {
        {[1 2 3]},
        {[1 2 3],[1 2]},
        {[1 2 3],[1 2 3]},
        {[1 2 3],[1 2 3], [1]},
        {[1 2 3],[1 2 3], [1 3 4]},
    };

    out = [3,6,9,9,27];

    for i = 1:length(in)
        g = grid.Cartesian(in{i}{:});
        testCase.verifyEqual(g.N(),out(i));
    end
end


function testD(testCase)
    in  = {
        {[1 2 3]},
        {[1 2 3],[1 2]},
        {[1 2 3],[1 2 3]},
        {[1 2 3],[1 2 3], [1]},
        {[1 2 3],[1 2 3], [1 3 4]},
    };

    out = [1,2,2,3,3];

    for i = 1:length(in)
        g = grid.Cartesian(in{i}{:});
        testCase.verifyEqual(g.D(),out(i));
    end
end

function testSize(testCase)
    in  = {
        {[1 2 3]},
        {[1 2 3],[1 2]},
        {[1 2 3],[1 2 3]},
        {[1 2 3],[1 2 3], [1]},
        {[1 2 3],[1 2 3], [1 3 4]},
    };

    out = {
        [3],
        [3 2],
        [3 3],
        [3 3 1],
        [3 3 3],
    };

    for i = 1:length(in)
        g = grid.Cartesian(in{i}{:});
        testCase.verifyEqual(g.size(),out{i});
    end
end

function testPoints(testCase)
    in  = {
        {[1 2]},
        {[1 2],[3 4]},
        {[1 2],[3 4], [5 6]},
    };

    out = {
        [[1; 2]],
        [[1; 1; 2; 2],[3; 4; 3; 4]],
        [[1; 1; 1; 1; 2; 2; 2; 2],[3; 3; 4; 4; 3; 3; 4; 4],[ 5; 6; 5; 6; 5; 6; 5; 6]],
    };

    for i = 1:length(in)
        g = grid.Cartesian(in{i}{:});
        testCase.verifyEqual(g.points(),out{i});
    end
end

function testMatrices(testCase)
    in  = {
        {[1 2]},
        {[1 2],[3 4]},
        {[1 2],[3 4], [5 6]},
    };

    out{1}{1} = [1; 2];

    out{2}{1} = [1, 1; 2, 2];
    out{2}{2} = [3, 4; 3, 4];

    out{3}{1}(:,:,1) = [1, 1; 2, 2];
    out{3}{1}(:,:,2) = [1, 1; 2, 2];

    out{3}{2}(:,:,1) = [3, 4; 3, 4];
    out{3}{2}(:,:,2) = [3, 4; 3, 4];

    out{3}{3}(:,:,1) = [5, 5; 5, 5];
    out{3}{3}(:,:,2) = [6, 6; 6, 6];

    for i = 1:length(in)
        g = grid.Cartesian(in{i}{:});
        testCase.verifyEqual(g.matrices(),out{i});
    end
end


function testRestrictFuncInvalidInput(testCase)
    inG1  = {
        {[1 2 3 4 5]},
        {[1 2 3],[4 5 6 7 8]},
        {[1 2 3],[4 5 6 7 8]},
        {[1 2 3],[4 5 6 7 8]},
    };

    inG2  = {
        {[1 3 4 5]},
        {[1 3],[4 5 6 8]},
        {[1 3],[4 6 8]},
        {[1 3],[4 6 8]},
    };

    inGf = {
        [1; 2; 3; 4; 5],
        [14; 15; 16; 17; 18; 24; 25; 26; 27; 28; 34; 35; 36; 37; 38];
        [14; 15; 16; 17; 18; 24; 25; 26; 27; 28; 34; 35; 36];
        [14; 15; 16; 17; 18; 24; 25; 26; 27; 28; 34; 35; 36; 37; 38; 39; 40];
    };

    out = {
        'grid:Cartesian:restrictFunc:NonMatchingGrids',
        'grid:Cartesian:restrictFunc:NonMatchingGrids',
        'grid:Cartesian:restrictFunc:NonMatchingFunctionSize',
        'grid:Cartesian:restrictFunc:NonMatchingFunctionSize',
    };

    for i = 1:length(inG1)
        g1 = grid.Cartesian(inG1{i}{:});
        g2 = grid.Cartesian(inG2{i}{:});
        testCase.verifyError(@()g1.restrictFunc(inGf{i},g2),out{i});
    end
end

function testRestrictFunc(testCase)
    inG1  = {
        {[1 2 3 4 5]},
        {[1 2 3],[4 5 6 7 8]},
    };

    inG2  = {
        {[1 3 5]},
        {[1 3],[4 6 8]},
    };

    inGf = {
        [1; 2; 3; 4; 5],
        [14; 15; 16; 17; 18; 24; 25; 26; 27; 28; 34; 35; 36; 37; 38];
    };

    outGf = {
        [1; 3; 5],
        [14; 16; 18; 34; 36; 38];
    };

    for i = 1:length(inG1)
        g1 = grid.Cartesian(inG1{i}{:});
        g2 = grid.Cartesian(inG2{i}{:});
        testCase.verifyEqual(g1.restrictFunc(inGf{i}, g2), outGf{i});
    end
end

function testScaling(testCase)
    in = {[1 2 3], [1 2]};
    g = grid.Cartesian(in{:});

    testCase.verifyError(@()g.scaling(),'grid:Cartesian:NoScalingSet');

    g.h = [2 1];
    testCase.verifyEqual(g.scaling(),[2 1]);

end
