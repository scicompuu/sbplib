classdef Animation < handle
    properties
        timeStepper
        representationMaker
        updaters
    end

    % add input validation

    methods
        function obj = Animation(timeStepper, representationMaker, updaters);
            obj.timeStepper = timeStepper;
            obj.updaters = updaters;
            obj.representationMaker = representationMaker;
        end

        function update(obj, r)
            for i = 1:length(obj.updaters)
                obj.updaters{i}(r);
            end
            drawnow
        end

        function run(obj, tEnd, timeModifier, do_pause)
            default_arg('do_pause', false)

            function next_t = G(next_t)
                obj.timeStepper.evolve(next_t);
                r = obj.representationMaker(obj.timeStepper);
                obj.update(r);

                if do_pause
                    pause
                end
            end

            anim.animate(@G, obj.timeStepper.t, tEnd, timeModifier);
        end

        function step(obj, tEnd, do_pause)
            default_arg('do_pause', false)

            while obj.timeStepper.t < tEnd
                obj.timeStepper.step();

                r = obj.representationMaker(obj.timeStepper);
                obj.update(r);

                % TODO: Make it never go faster than a certain fram rate

                if do_pause
                    pause
                end
            end
        end

        function saveMovie(obj, tEnd, timeModifier, figureHandle, dirname)
            save_frame = anim.setup_fig_mov(figureHandle, dirname);

            function next_t = G(next_t)
                obj.timeStepper.evolve(next_t);
                r = obj.representationMaker(obj.timeStepper);
                obj.update(r);

                save_frame();
            end

            fprintf('Generating and saving frames to: ..\n')
            anim.animate(@G, obj.timeStepper.t, tEnd, timeModifier);
            fprintf('Generating movies...\n')
            cmd = sprintf('bash %s/+anim/make_movie.sh %s', sbplibLocation(),dirname);
            system(cmd);
        end
    end
end
