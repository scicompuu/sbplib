function done = timeTask(taskName)
    fprintf('%s', taskName);
    tStart = tic;

    function done_fun()
        fprintf(' - done %fs\n', toc(tStart));
    end
    done = @done_fun;
end