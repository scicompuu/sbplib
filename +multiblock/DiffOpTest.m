function tests = DiffOpTest()
    tests = functiontests(localfunctions);
end

function testCreation(testCase)
    do = newMultiblockOp();
end

function testSplitOp(testCase)
    testCase.verifyFail();
end

function testBoundary_condition(testCase)
    testCase.verifyFail();
end

function testInterface(testCase)
    testCase.verifyFail();
end

function testSize(testCase)
    mbDo = newMultiblockOp();
    testCase.verifyEqual(mbDo.size(), 15)
end


function do = mockDiffOp(size, bc, interface)
    do.H = 1;
    do.D = 1;

    do.size = size;
    do.boundary_condition = bc;
    do.interface = interface;
end


function do = newMultiblockOp()
    grids = {
        grid.Cartesian([0 1 2], [3 4 5]);
        grid.Cartesian([1 2 3], [10 20]);
    };

    conn = cell(2,2);
    conn{1, 2} = {'s','n'};

    mbGrid = multiblock.Grid(grids, conn);

    function [c, p] = boundary_condition(~,~,~,~)
        c = 1; p = 1;
    end

    function [c, p] = interface(~,~,~,~)
        c = 1; p = 1;
    end

    doHand = @(grid,~)mockDiffOp(@(~)prod(grid.size()), @boundary_condition, @interface);

    do = multiblock.DiffOp(doHand, mbGrid, 0);
end