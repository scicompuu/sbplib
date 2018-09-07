function tests = equidistantTest()
    tests = functiontests(localfunctions);
end


function testErrorInvalidLimits(testCase)
     in  = {
        {10,{1}},
        {10,[0,1]},
        {[10,10],{0,1},{1}},
        {[10,10],{1},{1,0}},
        {10,{1,0}},
        {[10, 5],{1,0}, {0,-1}},
    };

    for i = 1:length(in)
        testCase.verifyError(@()grid.equidistant(in{i}{:}),'grid:equidistant:InvalidLimits',sprintf('in(%d) = %s',i,toString(in{i})));
    end
end

function testErrorNonMatchingParam(testCase)
    in  = {
        {[],{1}},
        {[],{1},{0,1}},
        {[5,5],{0,1},{0,1},{0,1}},
        {[5,5,4],{0,1},{0,1}},
    };

    for i = 1:length(in)
        testCase.verifyError(@()grid.equidistant(in{i}{:}),'grid:equidistant:NonMatchingParameters',sprintf('in(%d) = %s',i,toString(in{i})));
    end
end


function testCompiles(testCase)
    in  = {
        {5, {0,1}},
        {[3 3],{0,1},{0,2}},
    };

    out = {
        [[0; 0.25; 0.5; 0.75; 1]],
        [[0; 0; 0; 0.5; 0.5; 0.5; 1; 1; 1;],[0; 1; 2; 0; 1; 2; 0; 1; 2;]],
    };

    for i = 1:length(in)
        g = grid.equidistant(in{i}{:});
        testCase.verifyEqual(g.points(),out{i});
    end
end
