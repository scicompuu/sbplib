classdef Rungekutta4proper < time.Timestepper
    properties
        F
        k
        t
        v
        m
        n
    end


    methods
        % Timesteps v_t = F(v,t), using RK4 fromt t = t0 with timestep k and initial conditions v = v0
        function obj = Rungekutta4proper(F, k, t0, v0)
            obj.F = F;
            obj.k = k;
            obj.t = t0;
            obj.v = v0;
            obj.m = length(v0);
            obj.n = 0;
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