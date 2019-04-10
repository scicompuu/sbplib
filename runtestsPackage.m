function res = runtestsPackage(pkgName)
    ts = matlab.unittest.TestSuite.fromPackage(pkgName, 'IncludingSubpackages', true);
    res = ts.run();
end
