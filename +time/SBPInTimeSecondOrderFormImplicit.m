classdef SBPInTimeSecondOrderFormImplicit < time.Timestepper
    properties
        A, B, C, f
        AA, BB, ff

        n
        t
        k

        firstOrderTimeStepper
    end

    methods
        % Solves A*u_tt + B*u_t + C*u = f(t)
        % A, B can either both be constants or both be function handles,
        % They can also be omitted by setting them equal to the empty matrix.
        function obj = SBPInTimeSecondOrderFormImplicit(A, B, C, f, k, t0, v0, v0t, do_scaling, TYPE, order, blockSize)
            default_arg('f', []);
            default_arg('TYPE', []);
            default_arg('order', []);
            default_arg('blockSize',[]);
            default_arg('do_scaling', true);

            m = length(v0);

            default_arg('A', sparse(m, m));
            default_arg('B', sparse(m, m));
            default_arg('C', sparse(m, m));

            I = speye(m);
            O = sparse(m,m);

            % Rewrite to
            % AA*w_t = BB*w + ff(t);

            obj.AA = [
                 I, O;
                 O, A;
            ];
            obj.BB = [
                 O,  I;
                -C, -B;
            ];

            if ~isempty(f)
                obj.ff = @(t)[
                    sparse(m,1);
                           f(t);
                ];
            else
                obj.ff = @(t) sparse(2*m,1);
            end

            w0 = [v0; v0t];

            obj.k = k;
            obj.t = t0;
            obj.n = 0;

            if do_scaling
                scaling = [ones(m,1); sqrt(diag(C))];
                obj.firstOrderTimeStepper = time.SBPInTimeScaled(obj.AA, obj.BB, obj.ff, obj.k, obj.t, w0, scaling, TYPE, order, blockSize);
            else
                obj.firstOrderTimeStepper = time.SBPInTimeImplicitFormulation(obj.AA, obj.BB, obj.ff, obj.k, obj.t, w0, TYPE, order, blockSize);
            end
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
