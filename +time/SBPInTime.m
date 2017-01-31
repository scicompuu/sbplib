classdef SBPInTime < time.Timestepper
    % The SBP in time method.
    % Implemented for v_t = A*v + f(t)
    % k_local -- time-step
    % Nblock -- number of points in each block
    % nodes -- points such that t_n + nodes are the points in block n.
    % Each "step" takes one block step and thus advances
    % k = k_local*(Nblock-1) in time.
    % M -- matrix used in every solve.
    % [L,U,P,Q] = lu(M);
    properties
        M     % System matrix
        L,U,P % LU factorization of M
        Q
        A
        Et_r
        penalty
        f
        k_local
        k
        t
        v
        m
        n
        Nblock
        order
        nodes
    end

    methods
        function obj = SBPInTime(A, f, k, order, Nblock, t0, v0, TYPE)
            default_arg('TYPE','equidistant');
            default_arg('Nblock',time.SBPInTime.smallestBlockSize(order,TYPE));

            obj.A = A;
            obj.f = f;
            obj.k_local = k;
            obj.k = k*(Nblock-1);
            obj.Nblock = Nblock;
            obj.t = t0;
            obj.m = length(v0);
            obj.n = 0;

            %==== Build the time discretization matrix =====%
            switch TYPE
                case 'equidistant'
                    ops = sbp.D2Standard(Nblock,{0,obj.k},order);
                case 'optimal'
                    ops = sbp.D1Nonequidistant(Nblock,{0,obj.k},order);
                case 'minimal'
                    ops = sbp.D1Nonequidistant(Nblock,{0,obj.k},order,'minimal');
            end

            D1 = ops.D1;
            HI = ops.HI;
            e_l = ops.e_l;
            e_r = ops.e_r;
            obj.nodes = ops.x;

            Ix = speye(size(A));
            It = speye(Nblock,Nblock);

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
                              obj.penalty, obj.f, obj.Nblock,...
                              obj.Et_r,...
                              obj.L, obj.U, obj.P, obj.Q);
            obj.t = obj.t + obj.k;
            obj.n = obj.n + obj.Nblock-1;
        end
    end


    methods(Static)

        %
        function [k,numberOfBlocks] = alignedTimeStep(k,Tend,Nblock)

            % input k is the desired time-step
            % Nblock is the number of points per block.

            % Make sure that we reach the final time by advancing
            % an integer number of blocks
            kblock = (Nblock-1)*k;
            numberOfBlocks = ceil(Tend/kblock);
            kblock = Tend/(numberOfBlocks);

            % Corrected time step
            k = kblock/(Nblock-1);

        end

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
