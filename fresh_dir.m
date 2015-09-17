function fresh_dir(dirname)
    if isdir(dirname)
            fprintf('Directory %s already exists.\n',dirname);
            yN = input('Delete its content and continue? [y/N]: ','s');
            if strcmp(yN,'y')
                system(['rm -rf ' dirname]);
            else
                error('Can''t use directory.')
            end
        end
    mkdir(dirname);
end