classdef Timestepper < handle
    properties (Abstract)
        t
        k
        n
    end

    methods (Abstract)
         [v,t] = getV(obj)
         obj = step(obj)
    end


    methods
        function [v,t] = stepN(obj,n,progress_bar)
            default_arg('progress_bar',false);

            if ~progress_bar
                obj.stepN_without_progress(n);
            else
                obj.stepN_with_progress(n);
            end

            v = obj.getV;
            t = obj.t;
        end

        function stepN_without_progress(obj, n)
            for i=1:n
                obj.step();
            end
        end

        function stepN_with_progress(obj, n)
            FRAME_RATE = 20;

            steps_to_update = 1;
            steps_since_update = 0;
            last_update = tic();
            s = util.replace_string('','   %d %%',0);
            for i=1:n
                obj.step();

                steps_since_update = steps_since_update + 1;

                if steps_since_update >= steps_to_update
                    s = util.replace_string(s,'   %.2f %%',i/n*100);

                    time_since_update = toc(last_update);
                    time_error = time_since_update - 1/FRAME_RATE;
                    time_per_step = time_since_update/steps_since_update;

                    steps_to_update = max(steps_to_update - 0.9*time_error/time_per_step ,1);

                    steps_since_update = 0;
                    last_update = tic();
                end
            end

            s = util.replace_string(s,'');
        end

        function [v,t] = evolve(obj, tend, progress_bar)
            default_arg('progress_bar',false)
            if ~progress_bar
                obj.evolve_without_progress(tend);
            else
                obj.evolve_with_progress(tend);
            end
            v = obj.getV;
            t = obj.t;
        end

        function evolve_without_progress(obj, tend)
            while obj.t < tend - obj.k/100
                obj.step();
            end
        end

        function evolve_with_progress(obj, tend)
            FRAME_RATE = 20;

            steps_to_update = 1;
            steps_since_update = 0;
            last_update = tic();
            s = util.replace_string('','   %d %%',0);
            while obj.t < tend - obj.k/100
                obj.step();

                steps_since_update = steps_since_update + 1;

                if steps_since_update >= steps_to_update
                    s = util.replace_string(s,'   %.2f %%',obj.t/tend*100);

                    time_since_update = toc(last_update);
                    time_error = time_since_update - 1/FRAME_RATE;
                    time_per_step = time_since_update/steps_since_update;

                    steps_to_update = max(steps_to_update - 0.9*time_error/time_per_step ,1);

                    steps_since_update = 0;
                    last_update = tic();
                end
            end

            s = util.replace_string(s,'');
        end


    end
end