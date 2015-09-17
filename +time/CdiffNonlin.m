classdef CdiffNonlin < time.Timestepper
    properties
        D
        E
        S
        k
        t
        v
        v_prev
        n
    end


    methods
        function obj = CdiffNonlin(D, E, S, k, t0, v, v_prev)
            default_arg('S',0);
            default_arg('E',0);

            if isnumeric(S) && S == 0
                S = @(v)0;
            end

            if isnumeric(E) && E == 0
                E = @(v)0;
            end


            % m = size(D,1);
            % default_arg('E',sparse(m,m));
            % default_arg('S',sparse(m,1));

            obj.D = D;
            obj.E = E;
            obj.S = S;
            obj.k = k;
            obj.t = t0;
            obj.v = v;
            obj.v_prev = v_prev;
        end

        function [v,t] = getV(obj)
            v = obj.v;
            t = obj.t;
        end

        function [vt,t] = getVt(obj)
            vt = (obj.v-obj.v_prev)/obj.k; % Could be improved using u_tt = f(u))
            t = obj.t;
        end

        function obj = step(obj)
            [obj.v, obj.v_prev] = time.cdiff.cdiff(obj.v, obj.v_prev, obj.k, obj.D(obj.v), obj.E(obj.v), obj.S(obj.v));
            obj.t = obj.t + obj.k;
            obj.n = obj.n + 1;
        end
    end
end