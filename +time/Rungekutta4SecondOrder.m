classdef Rungekutta4SecondOrder < time.Timestepper
    properties
        F
        k
        t
        w
        m
        D
        E
        S
        M
        C
        n
    end


    methods
        % Solves u_tt = Du + Eu_t + S
        % D, E, S can either all be constants or all be function handles,
        % They can also be omitted by setting them equal to the empty matrix.
        function obj = Rungekutta4SecondOrder(D, E, S, k, t0, v0, v0t)
            obj.D = D;
            obj.E = E;
            obj.S = S;
            obj.m = length(v0);
            obj.n = 0;


            if isa(D, 'function_handle') || isa(E, 'function_handle') || isa(S, 'function_handle')
                default_arg('D', @(t)sparse(obj.m, obj.m));
                default_arg('E', @(t)sparse(obj.m, obj.m));
                default_arg('S', @(t)sparse(obj.m, 1)    );

                if ~isa(D, 'function_handle')
                    D = @(t)D;
                end
                if ~isa(D, 'function_handle')
                    E = @(t)E;
                end
                if ~isa(D, 'function_handle')
                    S = @(t)S;
                end

                I = speye(obj.m);
                O = sparse(obj.m,obj.m);

                obj.M = @(t)[
                       O,    I;
                    D(t), E(t);
                ];
                obj.C = @(t)[
                    zeros(obj.m,1);
                              S(t);
                ];

                obj.k = k;
                obj.t = t0;
                obj.w = [v0; v0t];

                obj.F = @(w,t)(obj.M(t)*w + obj.C(t));
            else

                default_arg('D', sparse(obj.m, obj.m));
                default_arg('E', sparse(obj.m, obj.m));
                default_arg('S', sparse(obj.m, 1)    );

                I = speye(obj.m);
                O = sparse(obj.m,obj.m);

                obj.M = [
                    O, I;
                    D, E;
                ];
                obj.C = [
                    zeros(obj.m,1);
                                 S;
                ];

                obj.k = k;
                obj.t = t0;
                obj.w = [v0; v0t];

                obj.F = @(w,t)(obj.M*w + obj.C);
            end
        end

        function [v,t] = getV(obj)
            v = obj.w(1:end/2);
            t = obj.t;
        end

        function [vt,t] = getVt(obj)
            vt = obj.w(end/2+1:end);
            t = obj.t;
        end

        function obj = step(obj)
            obj.w = time.rk4.rungekutta_4(obj.w, obj.t, obj.k, obj.F);
            obj.t = obj.t + obj.k;
            obj.n = obj.n + 1;
        end
    end


    methods (Static)
        function k = getTimeStep(lambda)
            k = rk4.get_rk4_time_step(lambda);
        end
    end

end