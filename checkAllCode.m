% Run matlabs checkcode on all files in current folder including subfolders
function checkAllCode(d)
    default_arg('d','')
    files = findMfiles(d);
    checkcode(files)
end
