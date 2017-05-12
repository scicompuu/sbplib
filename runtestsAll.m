function res = runtestsAll()
    l = dir();

    warning('Simplify using the ''what'' command')

    dirNames = {l([l.isdir]).name};

    packages = {};
    for i = 1:length(dirNames)
        if dirNames{i}(1) == '+'
            packages{end+1} = dirNames{i}(2:end);
        end
    end

    rootSuite = matlab.unittest.TestSuite.fromFolder(pwd);
    packageSuites = {};
    for i = 1:length(packages)
        packageSuites{i} = matlab.unittest.TestSuite.fromPackage(packages{i});
    end

    ts = [rootSuite, packageSuites{:}];

    res = ts.run();
end
