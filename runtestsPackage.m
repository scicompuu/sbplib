function res = runtestsPackage(pkgName)
    ts = matlab.unittest.TestSuite.fromPackage(pkgName);
    res = ts.run();
end