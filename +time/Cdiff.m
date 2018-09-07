classdef Cdiff < time.Timestepper
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
        % Solves u_tt = Du + Eu_t + S
        % D, E, S can either all be constants or all be function handles,
        % They can also be omitted by setting them equal to the empty matrix.
        % Cdiff(D, E, S, k, t0, n0, v, v_prev)
        function obj = Cdiff(D, E, S, k, t0, n0, v, v_prev)
            m = length(v);
            default_arg('E',sparse(m,m));
            default_arg('S',sparse(m,1));

            obj.D = D;
            obj.E = E;
            obj.S = S;


            obj.k = k;
            obj.t = t0;
            obj.n = n0;
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
            [obj.v, obj.v_prev] = time.cdiff.cdiff(obj.v, obj.v_prev, obj.k, obj.D, obj.E, obj.S);
            obj.t = obj.t + obj.k;
            obj.n = obj.n + 1;
        end
    end
end