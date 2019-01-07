classdef Cdiff < time.Timestepper
    properties
        A, B, C
        AA, BB, CC
        G
        k
        t
        v
        v_prev
        n
    end


    methods
        % Solves
        %   A*v_tt + B*v_t + C*v = G(t)
        %   v(t0) = v0
        %   v_t(t0) = v0t
        % starting at time t0 with timestep k
        % Using
        % A*Dp*Dm*v_n + B*D0*v_n + C*v_n = G(t_n)
        function obj = Cdiff(A, B, C, G, v0, v0t, k, t0)
            m = length(v0);
            default_arg('A', speye(m));
            default_arg('B', sparse(m,m));
            default_arg('G', @(t) sparse(m,1));
            default_arg('t0', 0);

            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.G = G;

            % Rewrite as AA*v_(n+1) + BB*v_n + CC*v_(n-1) = G(t_n)
            obj.AA = A/k^2 + B/2/k;
            obj.BB = -2*A/k^2 + C;
            obj.CC = A/k^2 - B/2/k;

            obj.k = k;
            obj.v_prev = v0;
            obj.v = v0 + k*v0t;
            obj.t = t0+k;
            obj.n = 1;
        end

        function [v,t] = getV(obj)
            v = obj.v;
            t = obj.t;
        end

        function [vt,t] = getVt(obj)
            vt = (obj.v-obj.v_prev)/obj.k; % Could be improved using u_tt = f(u))
            t = obj.t;
        end

        function E = getEnergy(obj)
            v  = obj.v;
            vp = obj.v_prev;
            vt = (obj.v - obj.v_prev)/obj.k;

            E = vt'*obj.A*vt + v'*obj.C*vp;
        end

        function obj = step(obj)
            v_next = obj.AA\(-obj.BB*obj.v - obj.CC*obj.v_prev + obj.G(obj.t));

            obj.v_prev = obj.v;
            obj.v      = v_next;
            obj.t = obj.t + obj.k;
            obj.n = obj.n + 1;
        end
    end
end