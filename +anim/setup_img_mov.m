function save_frame = setup_img_mov(dirname,filename_prefix,ulim,map,F)
    default_arg('map',parula(256));
    default_arg('F',@(x)x);

    fresh_dir(dirname);

    n = size(map,1);

    frame_nr = 0;
    function save_frame_F(Z)
        filename = sprintf('%s/%s%04d.png',dirname,filename_prefix, frame_nr);
        z = uint8(F(Z)/ulim(2)*n);
        imwrite(z,map,filename);
        frame_nr = frame_nr + 1;
    end

    save_frame = @save_frame_F;
end
