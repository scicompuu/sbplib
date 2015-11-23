classdef Rungekutta4 < time.Timestepper
    properties
        D
        S
        F
        k
        t
        v
        m
        n
    end


    methods
        function obj = Rungekutta4(D, S, k, t0, v0)
            obj.D = D;
            obj.k = k;
            obj.t = t0;
            obj.v = v0;
            obj.m = length(v0);
            obj.n = 0;

            if S == 0
                obj.S = zeros(obj.m,1);
            else
                obj.S = S;
            end

            if S == 0
                obj.F = @(v,t)(obj.D*v);
            else
                obj.F = @(v,t)(obj.D*v + obj.S);
            end
        end

        function [v,t] = getV(obj)
            v = obj.v;
            t = obj.t;
        end

        function obj = step(obj)
            obj.v = time.rk4.rungekutta_4(obj.v, obj.t, obj.k, obj.F);
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