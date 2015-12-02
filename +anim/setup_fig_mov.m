function save_frame = setup_fig_mov(figure_handle,dirname)

    fresh_dir(dirname)

    figure_handle.Units = 'centimeters';
    figure_handle.PaperUnits = 'centimeters';
    figure_handle.PaperPosition(3:4) = figure_handle.Position(3:4);


    frame_nr = 0;
    function save_frame_F()
        saveas(figure_handle,sprintf('%s/%s%04d.png',dirname,'fig',frame_nr));
        frame_nr = frame_nr + 1;
    end
    save_frame = @save_frame_F;
end
