function tests = DiffOpTest()
    tests = functiontests(localfunctions);
end

function testCreation(testCase)
    g = multiblock.Grid({},{});
    doHand = @(grid,order)[];
    order = 0;
    do = multiblock.DiffOp(doHand, g, order);
end




% function do = mockDiffOp()
%     do.H = 1;
%     do.D = 1;
% end