function [files,res] = checkSbplib()


    files = collectTargets([]);

    if nargout == 0
        checkcode(files, '-config=checksettings.txt');
        return
    end

    res = checkcode(files, '-config=checksettings.txt');

    % Remove any empty respones
    I = [];
    for i = 1:length(res)
        if isempty(res{i})
            I(end+1) = i;
        end
    end

    files(I) = [];
    res(I) = [];
end

function targets = collectTargets(dirPath)
    [mfiles, packages] = getFilesAndPackages(dirPath);

    targets = {};
    for i = 1:length(mfiles)
        targets{i} = fullfile(dirPath, mfiles{i});
    end

    for i = 1:length(packages)
        subtargets = collectTargets(fullfile(dirPath, packages{i}));
        targets = [targets subtargets];
    end
end

function [mfiles, packages] = getFilesAndPackages(dirPath)
    if isempty(dirPath)
        l = dir();
    else
        l = dir(dirPath);
    end

    packages = {};
    mfiles = {};

    for i = 1:length(l)
        if l(i).isdir && l(i).name(1) == '+'
            packages{end+1} = l(i).name;
        elseif ~l(i).isdir && strcmp(l(i).name(end-1:end),'.m')
            mfiles{end+1} = l(i).name;
        end
    end
end