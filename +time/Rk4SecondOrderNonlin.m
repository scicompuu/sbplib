classdef Rk4SecondOrderNonlin < time.Timestepper
    properties
        F
        k
        t
        w
        m

        D
        E
        S

        n
    end


    methods
        function obj = Rk4SecondOrderNonlin(D, E, S, k, t0, v0, v0t)

            if S == 0
                S = @(v,t)0;
            end

            if E == 0
                E = @(v,t)0;
            end

            obj.k = k;
            obj.t = t0;
            obj.w = [v0; v0t];

            m = length(v0);
            function wt = F(w,t)
                v  = w(1:m);
                vt = w(m+1:end);

                % Def: w = [v; vt]
                wt(1:m,1) = vt;
                wt(m+1:2*m,1) = D(v)*v + E(v)*vt + S(v,t);

            end

            obj.F = @F;
            obj.D = D;
            obj.E = E;
            obj.S = S;
            obj.m = m;
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