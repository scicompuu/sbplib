classdef SBPInTimeImplicitFormulation < time.Timestepper
    % The SBP in time method.
    % Implemented for A*v_t = B*v + f(t), v(0) = v0
    properties
        A,B
        f

        k % total time step.

        blockSize % number of points in each block
        N % Number of components

        order
        nodes

        M,K     % System matrices
        L,U,p,q % LU factorization of M
        e_T

        % Time state
        t
        v
        n
    end

    methods
        function obj = SBPInTimeImplicitFormulation(A, B, f, k, t0, v0, TYPE, order, blockSize)

            default_arg('TYPE','gauss');

            if(strcmp(TYPE,'gauss'))
                default_arg('order',4)
                default_arg('blockSize',4)
            else
                default_arg('order', 8);
                default_arg('blockSize',time.SBPInTimeImplicitFormulation.smallestBlockSize(order,TYPE));
            end

            obj.A = A;
            obj.B = B;
            obj.f = f;

            obj.k = k;
            obj.blockSize = blockSize;
            obj.N = length(v0);

            obj.n = 0;
            obj.t = t0;

            %==== Build the time discretization matrix =====%
            switch TYPE
                case 'equidistant'
                    ops = sbp.D2Standard(blockSize,{0,obj.k},order);
                case 'optimal'
                    ops = sbp.D1Nonequidistant(blockSize,{0,obj.k},order);
                case 'minimal'
                    ops = sbp.D1Nonequidistant(blockSize,{0,obj.k},order,'minimal');
                case 'gauss'
                    ops = sbp.D1Gauss(blockSize,{0,obj.k});
            end

            I = speye(size(A));
            I_t = speye(blockSize,blockSize);

            D1 = kron(ops.D1, I);
            HI = kron(ops.HI, I);
            e_0 = kron(ops.e_l, I);
            e_T = kron(ops.e_r, I);
            obj.nodes = ops.x;

            % Convert to form M*w = K*v0 + f(t)
            tau = kron(I_t, A) * e_0;
            M = kron(I_t, A)*D1 + HI*tau*e_0' - kron(I_t, B);

            K = HI*tau;

            obj.M = M;
            obj.K = K;
            obj.e_T = e_T;

            % LU factorization
            [obj.L,obj.U,obj.p,obj.q] = lu(obj.M, 'vector');

            obj.v = v0;
        end

        function [v,t] = getV(obj)
            v = obj.v;
            t = obj.t;
        end

        function obj = step(obj)
            RHS = zeros(obj.blockSize*obj.N,1);

            for i = 1:length(obj.blockSize)
                RHS((1 + (i-1)*obj.N):(i*obj.N)) = obj.f(obj.nodes(i));
            end

            RHS = RHS + obj.K*obj.v;

            y = obj.L\RHS(obj.p);
            z = obj.U\y;

            w = zeros(size(z));
            w(obj.q) = z;

            obj.v = obj.e_T'*w;

            obj.t = obj.t + obj.k;
            obj.n = obj.n + 1;
        end
    end

    methods(Static)
        function N = smallestBlockSize(order,TYPE)
            default_arg('TYPE','gauss')

            switch TYPE
                case 'gauss'
                    N = 4;
            end
        end
    end
end
