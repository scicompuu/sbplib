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

            if progress_bar && n > 1500
                n1000 = floor(n/1000);

                s = util.replace_string('','   %d %%',0);
                for i=1:n
                    obj.step();
                    if mod(i,n1000) == 0
                    s = util.replace_string(s,'   %.2f %%',i/n*100);
                    end
                end
                s = util.replace_string(s,'');
            else
                for i=1:n
                    obj.step();
                end
            end
            v = obj.getV;
            t = obj.t;
        end

        function [v,t] = evolve(obj, tend, progress_bar)
            default_arg('progress_bar',false)
            if progress_bar
                obj.evolve_with_progress(tend);
            else
                obj.evolve_without_progress(tend);
            end
            v = obj.getV;
            t = obj.t;
        end

        function evolve_with_progress(obj, tend)
            dt = tend-obj.t;
            n = floor(dt/obj.k);
            n1000 = floor(n/1000);

            i = 0;
            s = util.replace_string('','   %d %%',0);
            while obj.t < tend - obj.k/100
                obj.step();

                i = i + 1;
                if mod(i,n1000) == 0
                    s = util.replace_string(s,'   %.2f %%',i/n*100);
                end
            end

            % if t < tend
            %     v = rk4.rungekutta_4(v, t, tend-t,F);
            % end

            s = util.replace_string(s,'');
        end

        function evolve_without_progress(obj, tend)
            while obj.t < tend - obj.k/100
                obj.step();
            end

            % if t < tend
            %     v = rk4.rungekutta_4(v, t, tend-t,F);
            % end
        end
    end
end