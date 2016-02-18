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