function prof(f)
    profile on
    try
        f();
        profile viewer
    catch e
        fprintf(2, '\n%s', getReport(e));
        profile clear
    end
end