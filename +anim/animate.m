% Calls F(t) repeatedly
% Should there be a Fsetup and a F, two function, to allow creating a plot and then updating it?
% F takes the time to generate the frame for and returns the actual time for the generated frame.
% t = F(t_r) is a function that paints a frame for time t. t is the closest time <=t_r
% it will be called for increasnig t.

%Todo: make it catch up and produce warnings if it lags behind? Instead of just requesting the next target time
function  animate(F, tstart, tend, time_modifier , frame_rate)
    if ~exist('time_modifier')
        time_modifier = 1;
    end

    if ~exist('frame_rate')
        frame_rate = 30;
    end

    frame_time = 1/frame_rate;
    dt = frame_time*time_modifier;

    animation_start = tic();
    t = F(tstart);
    while t < tend
        t = F(t + dt);
        t_left = (t-tstart)/time_modifier-toc(animation_start);
        pause(t_left)
    end
    time_to_animate = toc(animation_start);
    expected_time = tend/time_modifier;

    fprintf('\n');
    fprintf('Time to animate: %.3f\n', time_to_animate)
    fprintf('Expected time  : %.3f\n', expected_time)
end
