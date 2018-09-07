function tests = reshapeRowMajTest()
    tests = functiontests(localfunctions);
end

function test1D(testCase)
    in = {
        {5,[1; 2; 3; 4; 5]},
        {5,[1 2 3 4 5]},
    };
    out = {
        [1; 2; 3; 4; 5],
        [1; 2; 3; 4; 5],
    };
    for i = 1:length(in)
        testCase.verifyEqual(reshapeRowMaj(in{i}{2}, in{i}{1}),out{i});
    end
end


function testIdentity(testCase)
     in = {
        {[2,2], magic(2)},
        {[3,3], magic(3)},
        {[2,3], [1 2 3; 4 5 6]},
    };

    for i = 1:length(in)
        testCase.verifyEqual(reshapeRowMaj(in{i}{2}, in{i}{1}),in{i}{2});
    end
end

function test2D(testCase)
    in = {
        {[2,2],[11; 12; 21; 22]},
        {[3,2],[1 2 3; 4 5 6]},
        {[6 1],[1 2 3; 4 5 6]},
        {[1 6],[1 2 3; 4 5 6]},
    };

    out{1}(1,1) = 11;
    out{1}(1,2) = 12;
    out{1}(2,1) = 21;
    out{1}(2,2) = 22;

    out{2} = [1 2; 3 4; 5 6];
    out{3} = [1; 2; 3; 4; 5; 6];
    out{4} = [1 2 3 4 5 6];

    for i = 1:length(in)
        testCase.verifyEqual(reshapeRowMaj(in{i}{2}, in{i}{1}),out{i});
    end
end

function test3D(testCase)
    in = {
        {[2, 2, 2], [111; 112; 121; 122; 211; 212; 221; 222]},
        {[8 1], cat(3,[1 2; 3 4],[5 6; 7 8])},
        {[1 8], cat(3,[1 2; 3 4],[5 6; 7 8])},
        {[2 4], cat(3,[1 2; 3 4],[5 6; 7 8])},
        {[4 2], cat(3,[1 2; 3 4],[5 6; 7 8])},
    };

    out{1}(1,1,1) = 111;
    out{1}(1,1,2) = 112;
    out{1}(1,2,1) = 121;
    out{1}(1,2,2) = 122;
    out{1}(2,1,1) = 211;
    out{1}(2,1,2) = 212;
    out{1}(2,2,1) = 221;
    out{1}(2,2,2) = 222;

    out{2} = [1; 5; 2; 6; 3; 7; 4; 8];
    out{3} = [1  5  2  6  3  7  4  8];
    out{4} = [1  5  2  6;  3  7  4  8];
    out{5} = [1  5;  2  6;  3  7;  4  8];

    for i = 1:length(in)
        testCase.verifyEqual(reshapeRowMaj(in{i}{2}, in{i}{1}),out{i});
    end
end

function testNonSquare(testCase)
    in = {
        {[2, 3, 4],[111; 112; 113; 114; 121; 122; 123; 124; 131; 132; 133; 134; 211; 212; 213; 214; 221; 222; 223; 224; 231; 232; 233; 234]},
    };
    out{1}(1,1,1) = 111;
    out{1}(1,1,2) = 112;
    out{1}(1,1,3) = 113;
    out{1}(1,1,4) = 114;
    out{1}(1,2,1) = 121;
    out{1}(1,2,2) = 122;
    out{1}(1,2,3) = 123;
    out{1}(1,2,4) = 124;
    out{1}(1,3,1) = 131;
    out{1}(1,3,2) = 132;
    out{1}(1,3,3) = 133;
    out{1}(1,3,4) = 134;
    out{1}(2,1,1) = 211;
    out{1}(2,1,2) = 212;
    out{1}(2,1,3) = 213;
    out{1}(2,1,4) = 214;
    out{1}(2,2,1) = 221;
    out{1}(2,2,2) = 222;
    out{1}(2,2,3) = 223;
    out{1}(2,2,4) = 224;
    out{1}(2,3,1) = 231;
    out{1}(2,3,2) = 232;
    out{1}(2,3,3) = 233;
    out{1}(2,3,4) = 234;
    for i = 1:length(in)
        testCase.verifyEqual(reshapeRowMaj(in{i}{2}, in{i}{1}),out{i});
    end
end