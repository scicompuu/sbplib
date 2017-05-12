% Find all m files reachable from folder d
function files = findMfiles(d)
    contents = what(d);

    files = {};
    for i = 1:length(contents.m)
        files{end+1} = fullfile(d, contents.m{i});
    end

    for i = 1:length(contents.packages)
        packagePath = fullfile(d,['+' contents.packages{i}]);
        files = [files, findMfiles(packagePath)];
    end
end