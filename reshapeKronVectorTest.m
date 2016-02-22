function tests = reshapeKronVectorTest()
    tests = functiontests(localfunctions);
end

function test1D(testCase)
    inGf = [1 2 3 4 5]';
    inM = 5;
    out = [1 2 3 4 5]';
    testCase.verifyEqual(reshapeKronVector(inGf, inM),out);
end

function test2D(testCase)
    inGf = [11; 12; 21; 22];
    inM = [2, 2];

    out(1,1) = 11;
    out(1,2) = 12;
    out(2,1) = 21;
    out(2,2) = 22;

    testCase.verifyEqual(reshapeKronVector(inGf, inM),out);
end

function test3D(testCase)
    inGf = [111; 112; 121; 122; 211; 212; 221; 222];
    inM = [2, 2, 2];

    out(1,1,1) = 111;
    out(1,1,2) = 112;
    out(1,2,1) = 121;
    out(1,2,2) = 122;
    out(2,1,1) = 211;
    out(2,1,2) = 212;
    out(2,2,1) = 221;
    out(2,2,2) = 222;

    testCase.verifyEqual(reshapeKronVector(inGf, inM),out);
end

function testNonSquare(testCase)
    inGf = [
        111;
        112;
        113;
        114;
        121;
        122;
        123;
        124;
        131;
        132;
        133;
        134;
        211;
        212;
        213;
        214;
        221;
        222;
        223;
        224;
        231;
        232;
        233;
        234;
    ];
    inM = [2, 3, 4];

    out(1,1,1) = 111;
    out(1,1,2) = 112;
    out(1,1,3) = 113;
    out(1,1,4) = 114;
    out(1,2,1) = 121;
    out(1,2,2) = 122;
    out(1,2,3) = 123;
    out(1,2,4) = 124;
    out(1,3,1) = 131;
    out(1,3,2) = 132;
    out(1,3,3) = 133;
    out(1,3,4) = 134;
    out(2,1,1) = 211;
    out(2,1,2) = 212;
    out(2,1,3) = 213;
    out(2,1,4) = 214;
    out(2,2,1) = 221;
    out(2,2,2) = 222;
    out(2,2,3) = 223;
    out(2,2,4) = 224;
    out(2,3,1) = 231;
    out(2,3,2) = 232;
    out(2,3,3) = 233;
    out(2,3,4) = 234;

    testCase.verifyEqual(reshapeKronVector(inGf, inM), out);
end