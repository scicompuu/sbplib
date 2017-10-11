function location = sbplibLocation()
    scriptname  = mfilename('fullpath');
    [location, ~, ~] = fileparts(scriptname);
end
