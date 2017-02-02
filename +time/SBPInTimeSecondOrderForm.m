classdef SBPInTimeSecondOrderForm < time.Timestepper
    properties
        A,B,C
        M, f

        n
        t
        k

        firstOrderTimeStepper
    end

    methods
        % Solves u_tt = Au + Bu_t + C
        % A, B can either both be constants or both be function handles,
        % They can also be omitted by setting them equal to the empty matrix.
        function obj = SBPInTimeSecondOrderForm(A, B, C, k, t0, v0, v0t, TYPE, order, blockSize)
            default_arg('TYPE', []);
            default_arg('order', []);
            default_arg('blockSize',[]);

            m = length(v0);

            default_arg('A', sparse(m, m));
            default_arg('B', sparse(m, m));
            default_arg('C', sparse(m, 1));

            I = speye(m);
            O = sparse(m,m);

            obj.M = [
                O, I;
                A, B;
            ];
            obj.f = @(t)[
                sparse(m,1);
                          C;
            ];

            w0 = [v0; v0t];

            obj.k = k;
            obj.t = t0;
            obj.n = 0;

            obj.firstOrderTimeStepper = time.SBPInTime(obj.M, obj.f, obj.k, obj.t, w0, TYPE, order, blockSize);
        end

        function [v,t] = getV(obj)
            w = obj.firstOrderTimeStepper.getV();
            v = w(1:end/2);
            t = obj.t;
        end

        function [vt,t] = getVt(obj)
            w = obj.firstOrderTimeStepper.getV();
            vt = w(end/2+1:end);
            t = obj.t;
        end

        function obj = step(obj)
            obj.firstOrderTimeStepper.step();
            obj.t = obj.firstOrderTimeStepper.t;
            obj.n = obj.firstOrderTimeStepper.n;
        end
    end
end
