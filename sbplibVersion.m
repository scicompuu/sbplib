% Prints the version and location of the sbplib currently in use.
function sbplibVersion()
    scriptname  = mfilename('fullpath');
    I = find(scriptname == '/',1,'last');
    folder = scriptname(1:I);

    name = 'sbplib';
    ver = '0.0.x';

    fprintf('%s %s\n', name, ver);
    fprintf('Running in:\n%s\n',folder);
end