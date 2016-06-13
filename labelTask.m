function o = labelTask(in)

    switch class(in)
        case 'char'
            fprintf(in);
            o = tic();
        case 'uint64'
            o = toc(in);
            fprintf(' - done %fs\n', o);
        otherwise
            error('Unknow input type: %s', class(in))
    end
