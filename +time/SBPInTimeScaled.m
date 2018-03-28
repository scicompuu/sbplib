classdef SBPInTimeScaled < time.Timestepper
    % The SBP in time method.
    % Implemented for A*v_t = B*v + f(t), v(0) = v0
    % The resulting system of equations is
    %   M*u_next= K*u_prev_end + f
    properties
        A,B
        f

        k % total time step.

        blockSize % number of points in each block
        N % Number of components

        order
        nodes

        Mtilde,Ktilde     % System matrices
        L,U,p,q % LU factorization of M
        e_T

        scaling
        S, Sinv % Scaling matrices

        % Time state
        t
        vtilde
        n
    end

    methods
        function obj = SBPInTimeScaled(A, B, f, k, t0, v0, scaling, TYPE, order, blockSize)
            default_arg('TYPE','gauss');
            default_arg('f',[]);

            if(strcmp(TYPE,'gauss'))
                default_arg('order',4)
                default_arg('blockSize',4)
            else
                default_arg('order', 8);
                default_arg('blockSize',time.SBPInTimeImplicitFormulation.smallestBlockSize(order,TYPE));
            end

            obj.A = A;
            obj.B = B;
            obj.scaling = scaling;

            if ~isempty(f)
                obj.f = f;
            else
                obj.f = @(t)sparse(length(v0),1);
            end

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

            obj.S =    kron(I_t, spdiag(scaling));
            obj.Sinv = kron(I_t, spdiag(1./scaling));

            obj.Mtilde = obj.Sinv*M*obj.S;
            obj.Ktilde = obj.Sinv*K*spdiag(scaling);
            obj.e_T = e_T;


            % LU factorization
            [obj.L,obj.U,obj.p,obj.q] = lu(obj.Mtilde, 'vector');

            obj.vtilde = (1./obj.scaling).*v0;
        end

        function [v,t] = getV(obj)
            v = obj.scaling.*obj.vtilde;
            t = obj.t;
        end

        function obj = step(obj)
            forcing = zeros(obj.blockSize*obj.N,1);

            for i = 1:obj.blockSize
                forcing((1 + (i-1)*obj.N):(i*obj.N)) = obj.f(obj.t + obj.nodes(i));
            end

            RHS = obj.Sinv*forcing + obj.Ktilde*obj.vtilde;

            y = obj.L\RHS(obj.p);
            z = obj.U\y;

            w = zeros(size(z));
            w(obj.q) = z;

            obj.vtilde = obj.e_T'*w;

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
