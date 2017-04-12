% Prints the version and location of the sbplib currently in use.
function sbplibVersion()
    scriptname  = mfilename('fullpath');
    [folder,~,~] = fileparts(scriptname);

    name = 'sbplib (feature/grids)';
    ver = '0.0.x';

    fprintf('%s %s\n', name, ver);
    fprintf('Running in:\n%s\n',folder);
end