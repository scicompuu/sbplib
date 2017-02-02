classdef SBPInTime < time.Timestepper
    % The SBP in time method.
    % Implemented for v_t = A*v + f(t)
    %
    % Each "step" takes one block step and thus advances
    % k = k_local*(blockSize-1) in time.
    properties
        M     % System matrix
        L,U,P,Q % LU factorization of M
        A
        Et_r
        penalty
        f
        k_local % step size within a block
        k % Time size of a block  k/(blockSize-1) = k_local
        t
        v
        m
        n
        blockSize % number of points in each block
        order
        nodes
    end

    methods
        function obj = SBPInTime(A, f, k, t0, v0, TYPE, order, blockSize)
            default_arg('TYPE','minimal');
            default_arg('order', 8);
            default_arg('blockSize',time.SBPInTime.smallestBlockSize(order,TYPE));

            obj.A = A;
            obj.f = f;
            obj.k_local = k/(blockSize-1);
            obj.k = k;
            obj.blockSize = blockSize;
            obj.t = t0;
            obj.m = length(v0);
            obj.n = 0;

            %==== Build the time discretization matrix =====%
            switch TYPE
                case 'equidistant'
                    ops = sbp.D2Standard(blockSize,{0,obj.k},order);
                case 'optimal'
                    ops = sbp.D1Nonequidistant(blockSize,{0,obj.k},order);
                case 'minimal'
                    ops = sbp.D1Nonequidistant(blockSize,{0,obj.k},order,'minimal');
            end

            D1 = ops.D1;
            HI = ops.HI;
            e_l = ops.e_l;
            e_r = ops.e_r;
            obj.nodes = ops.x;

            Ix = speye(size(A));
            It = speye(blockSize,blockSize);

            obj.Et_r = kron(e_r,Ix);

            % Time derivative + penalty
            tau = 1;
            Mt = D1 + tau*HI*(e_l*e_l');

            % penalty to impose "data"
            penalty = tau*HI*e_l;
            obj.penalty = kron(penalty,Ix);

            Mx = kron(It,A);
            Mt = kron(Mt,Ix);
            obj.M = Mt - Mx;
            %==============================================%

            % LU factorization
            [obj.L,obj.U,obj.P,obj.Q] = lu(obj.M);

            % Pretend that the initial condition is the last level
            % of a previous step.
            obj.v = obj.Et_r * v0;

        end

        function [v,t] = getV(obj)
            v = obj.Et_r' * obj.v;
            t = obj.t;
        end

        function obj = step(obj)
            obj.v = time.sbp.sbpintime(obj.v, obj.t, obj.nodes,...
                              obj.penalty, obj.f, obj.blockSize,...
                              obj.Et_r,...
                              obj.L, obj.U, obj.P, obj.Q);
            obj.t = obj.t + obj.k;
            obj.n = obj.n + 1;
        end
    end


    methods(Static)
        function N = smallestBlockSize(order,TYPE)
            default_arg('TYPE','equidistant')

            switch TYPE

                case 'equidistant'
                    switch order
                        case 2
                            N = 2;
                        case 4
                            N = 8;
                        case 6
                            N = 12;
                        case 8
                            N = 16;
                        case 10
                            N = 20;
                        case 12
                            N = 24;
                        otherwise
                            error('Operator does not exist');
                    end

                case 'optimal'

                    switch order
                        case 4
                            N = 8;
                        case 6
                            N = 12;
                        case 8
                            N = 16;
                        case 10
                            N = 20;
                        case 12
                            N = 24;
                        otherwise
                            error('Operator does not exist');
                    end

                case 'minimal'

                    switch order
                        case 4
                            N = 6;
                        case 6
                            N = 10;
                        case 8
                            N = 12;
                        case 10
                            N = 16;
                        case 12
                            N = 20;
                        otherwise
                            error('Operator does not exist');
                    end
            end
        end
    end
end
