function tests = multiplyTest()
    tests = functiontests(localfunctions);
end


function testMultiply(testCase)
    a11 = [
        0.8147    0.1270;
        0.9058    0.9134;
    ];
    a12 = [
        0.6324    0.2785    0.9575;
        0.0975    0.5469    0.9649;
    ];
    a21 = [
        0.1576    0.9706;
    ];
    a22 = [
        0.9572    0.4854    0.8003;
    ];
    A = {
        a11 a12;
        a21 a22;
    };

    b11 = [
        0.1419    0.9157    0.9595;
        0.4218    0.7922    0.6557;
    ];
    b12 = [
        0.0357    0.9340;
        0.8491    0.6787;
    ];
    b21 = [
        0.7577    0.6555    0.0318;
        0.7431    0.1712    0.2769;
        0.3922    0.7060    0.0462;
    ];
    b22 = [
        0.0971    0.3171;
        0.8235    0.9502;
        0.6948    0.0344;
    ];

    B = {
        b11 b12;
        b21 b22;
    };


    C = {
        a11*b11 + a12*b21, a11*b12 + a12*b22;
        a21*b11 + a22*b21, a21*b12 + a22*b22;
    };

    testCase.verifyEqual(blockmatrix.multiply(A,B), C);
end