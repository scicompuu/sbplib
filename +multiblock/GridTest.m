function tests = GridTest()
    tests = functiontests(localfunctions);
end

function testCreation(testCase)
    g = multiblock.Grid({},{});
end

function testMissing(testCase)
    testCase.verifyFail();
end

function testGetBoundaryNames(testCase)
    [grids, conn] = prepareAdjecentBlocks();

    mbg = multiblock.Grid(grids, conn, multiblock.BoundaryGroup({1,'w'},{2,'w'}) );

    testCase.verifyFail();
end

function testGetBoundary(testCase)
    [grids, conn] = prepareAdjecentBlocks();

    mbg = multiblock.Grid(grids, conn, multiblock.BoundaryGroup({1,'w'},{2,'w'}) );
    testCase.verifyFail();
end


function [grids, conn] = prepareAdjecentBlocks()
    grids = {
        grid.Cartesian([0 1 2], [3 4 5]);
        grid.Cartesian([1 2], [10 20]);
    };

    conn = cell(2,2);
    conn{1, 2} = {'s','n'};
end