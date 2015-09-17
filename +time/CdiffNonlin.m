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
        function obj = Cdiff(D, E, S, k, t0, v, v_prev)
            m = size(D,1);
            default_arg('E',sparse(m,m));
            default_arg('S',sparse(m,1));

            if ~(issparse(D) && issparse(E) && issparse(S))
                warning('One of the matrices D, E, S is not sparse!')
                print_issparse(D)
                print_issparse(E)
                print_issparse(S)
            end

            obj.D = D;
            obj.E = E;
            obj.S = S;
            obj.k = k;
            obj.t = t0+k;
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