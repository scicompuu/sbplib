classdef CdiffImplicit < time.Timestepper
    properties
        A, B, C, G
        AA, BB, CC
        k
        t
        v, v_prev
        n

        % LU factorization
        L,U,p,q
    end

    methods
        % Solves
        %   A*u_tt + B*u + C*v_t = G(t)
        %   u(t0) = f1
        %   u_t(t0) = f2
        % starting at time t0 with timestep k
        function obj = CdiffImplicit(A, B, C, G, f1, f2, k, t0)
            default_arg('A', []);
            default_arg('C', []);
            default_arg('G', []);
            default_arg('f1', 0);
            default_arg('f2', 0);
            default_arg('t0', 0);

            m = size(B,1);

            if isempty(A)
                A = speye(m);
            end

            if isempty(C)
                C = sparse(m,m);
            end

            if isempty(G)
                G = @(t) sparse(m,1);
            end

            if isempty(f1)
                f1 = sparse(m,m);
            end

            if isempty(f2)
                f2 = sparse(m,m);
            end

            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.G = G;

            AA = 1/k^2*A + 1/2*B + 1/(2*k)*C;
            BB = -2/k^2*A;
            CC = 1/k^2*A + 1/2*B - 1/(2*k)*C;
            % AA*v_next + BB*v + CC*v_prev == G(t_n)

            obj.AA = AA;
            obj.BB = BB;
            obj.CC = CC;

            v_prev = f1;
            I = speye(m);
            v = (1/k^2*A)\((1/k^2*A - 1/2*B)*f1 + (1/k*I - 1/2*C)*f2 + 1/2*G(0));

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
            vt = (obj.v-obj.v_prev)/obj.k;
            t = obj.t;
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





%%% Derivation
% syms A B C G
% syms n k
% syms f1 f2

% v = symfun(sym('v(n)'),n);


% d = A/k^2 * (v(n+1) - 2*v(n) +v(n-1)) + B/2*(v(n+1)+v(n-1)) + C/(2*k)*(v(n+1) - v(n-1)) == G
% ic1 = v(0) == f1
% ic2 = A/k*(v(1)-f1) + k/2*(B*f1 + C*f2 - G) - f2 == 0

% c = collect(d, [v(n) v(n-1) v(n+1)]) % (-(2*A)/k^2)*v(n) + (B/2 + A/k^2 - C/(2*k))*v(n - 1) + (B/2 + A/k^2 + C/(2*k))*v(n + 1) == G
% syms AA BB CC
% % AA = B/2 + A/k^2 + C/(2*k)
% % BB = -(2*A)/k^2
% % CC = B/2 + A/k^2 - C/(2*k)
% s = subs(c, [B/2 + A/k^2 + C/(2*k), -(2*A)/k^2, B/2 + A/k^2 - C/(2*k)], [AA, BB, CC])


% ic2_a = collect(ic2, [v(1) f1 f2]) % (A/k)*v(1) + ((B*k)/2 - A/k)*f1 + ((C*k)/2 - 1)*f2 - (G*k)/2 == 0

