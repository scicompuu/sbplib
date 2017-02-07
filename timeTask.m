function done = timeTask(fmt, varargin)
    fprintf(fmt, varargin{:});
    tStart = tic;

    function done_fun()
        fprintf(' - done %fs\n', toc(tStart));
    end
    done = @done_fun;
end