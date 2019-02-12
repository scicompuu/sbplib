classdef CdiffImplicit < time.Timestepper
    properties
        A, B, C
        AA, BB, CC
        G
        k
        t
        v, v_prev
        n

        % LU factorization
        L,U,p,q
    end

    methods
        % Solves
        %   A*v_tt + B*v_t + C*v = G(t)
        %   v(t0) = v0
        %   v_t(t0) = v0t
        % starting at time t0 with timestep
        % Using
        % A*Dp*Dm*v_n + B*D0*v_n + C*I0*v_n = G(t_n)
        function obj = CdiffImplicit(A, B, C, G, v0, v0t, k, t0)
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
            AA =    A/k^2 + B/(2*k) + C/2;
            BB = -2*A/k^2;
            CC =    A/k^2 - B/(2*k) + C/2;

            obj.AA = AA;
            obj.BB = BB;
            obj.CC = CC;

            v_prev = v0;
            I = speye(m);
            v = v0 + k*v0t;

            if ~issparse(A) || ~issparse(B) || ~issparse(C)
                error('LU factorization with full pivoting only works for sparse matrices.')
            end

            [L,U,p,q] = lu(AA,'vector');

            obj.L = L;
            obj.U = U;
            obj.p = p;
            obj.q = q;

            obj.k = k;
            obj.t = t0+k;
            obj.n = 1;
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

        % Calculate the conserved energy (Dm*v_n)^2_A + Im*v_n^2_B
        function E = getEnergy(obj)
            v  = obj.v;
            vp = obj.v_prev;
            vt = (obj.v - obj.v_prev)/obj.k;

            E = vt'*obj.A*vt + 1/2*(v'*obj.C*v + vp'*obj.C*vp);
        end

        function obj = step(obj)
            b = obj.G(obj.t) - obj.BB*obj.v - obj.CC*obj.v_prev;
            obj.v_prev = obj.v;

            % % Backslash
            % obj.v = obj.AA\b;

            % LU with column pivot
            y = obj.L\b(obj.p);
            z = obj.U\y;
            obj.v(obj.q) = z;

            % Update time
            obj.t = obj.t + obj.k;
            obj.n = obj.n + 1;
        end
    end
end
