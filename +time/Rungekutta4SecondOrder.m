classdef Rungekutta4SecondOrder < time.Timestepper
    properties
        F
        k
        t
        w
        m
        D
        E
        S
        M
        C
        n
    end


    methods
        function obj = Rungekutta4SecondOrder(D, E, S, k, t0, v0, v0t)
            obj.D = D;
            obj.E = E;
            obj.S = S;
            obj.m = length(v0);

            I = speye(obj.m);
            O = sparse(obj.m,obj.m);
            obj.M = [O, I; D, E*I]; % Multiply with I to allow 0 as input.

            if S == 0
                obj.C = zeros(2*obj.m,1);
            else
                obj.C = [zeros(obj.m,1), S];
            end

            obj.k = k;
            obj.t = t0;
            obj.w = [v0; v0t];

            obj.F = @(w,t)(obj.M*w + obj.C);
        end

        function [v,t] = getV(obj)
            v = obj.w(1:end/2);
            t = obj.t;
        end

        function [vt,t] = getVt(obj)
            vt = obj.w(end/2+1:end);
            t = obj.t;
        end

        function obj = step(obj)
            obj.w = time.rk4.rungekutta_4(obj.w, obj.t, obj.k, obj.F);
            obj.t = obj.t + obj.k;
            obj.n = obj.n + 1;
        end
    end


    methods (Static)
        function k = getTimeStep(lambda)
            k = rk4.get_rk4_time_step(lambda);
        end
    end

end