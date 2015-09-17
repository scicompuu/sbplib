% Creates a unique filename on the form base_name.x.suffix where x is some number.
function filename = unique_filename(base_name, suffix)
    filename = strcat(base_name,suffix);
    i = 1;
    while exist(filename)
        filename = strcat(base_name,'.',num2str(i),suffix);
        i = i+1;
    end
end
