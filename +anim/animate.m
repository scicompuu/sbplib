% Calls F(t) repeatedly
% Should there be a Fsetup and a F, two function, to allow creating a plot and then updating it?
% F takes the time to generate the frame for and returns the actual time for the generated frame.
% t = F(t_r) is a function that paints a frame for time t. t is the closest time <=t_r
% it will be called for increasnig t.

%Todo: make it catch up and produce warnings if it lags behind? Instead of just requesting the next target time


% If adapt is true time_modifier is treated as an upper bound
function animate(F, tstart, tend, time_modifier, target_frame_rate)
    default_arg('time_modifier', 1);
    default_arg('target_frame_rate',30);

    % t is simulation time
    % tau is real time

    time_modifier_bound = time_modifier;
    dTau_target = 1/target_frame_rate; % Real time between frames

    rs = util.ReplaceableString();
    rs.appendFormat('              tau: %d\n');
    rs.appendFormat('       target tau: %d\n');
    rs.appendFormat('       Target fps: %.2f\n');
    rs.appendFormat('       Actual fps: %.2f\n');

    frameId = 0;
    NframeAvg = 20;
    frameTimes = ones(1,NframeAvg);

    animation_start = tic();
    prevTau = 0;
    targetTau = 0;
    t = F(tstart);

    while t < tend
        % Sleep until the frame should start
        pause(targetTau-toc(animation_start));
        tauFrameStart = targetTau;

        dt_target = dTau_target*time_modifier; % Targeted simulation time between frames

        t_prev = t;
        t = F(t + dt_target); % Run simulation


        % Calculate when this frame should end and the next start. (this depends on what simulation time we ended up on)
        dt = t-t_prev;
        targetTau = targetTau + dt/time_modifier;

        % Update information about this frame
        tau = toc(animation_start);
        rs.updateParam(tau, targetTau, 1/dTau_target, 1/(targetTau-tauFrameStart));
    end


    % Final time reporting
    time_to_animate = toc(animation_start);
    expected_time = tend/time_modifier;
    fprintf('\n');
    fprintf('Time to animate: %.3f\n', time_to_animate)
    fprintf('Expected time  : %.3f\n', expected_time)
end
