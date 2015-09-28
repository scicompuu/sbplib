classdef CdiffNonlin < time.Timestepper
    properties
        D
        E
        S
        k
        t
        v
        v_prev
        n
    end


    methods
        function obj = CdiffNonlin(D, E, S, k, t0,n0, v, v_prev)
            m = size(D(v),1);
            default_arg('E',0);
            default_arg('S',0);

            if isnumeric(S)
                S = @(v,t)S;
            end

            if isnumeric(E)
                E = @(v)E;
            end


            % m = size(D,1);
            % default_arg('E',sparse(m,m));
            % default_arg('S',sparse(m,1));

            obj.D = D;
            obj.E = E;
            obj.S = S;
            obj.k = k;
            obj.t = t0;
            obj.n = n0;
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

        function obj = step(obj)
            D = obj.D(obj.v);
            E = obj.E(obj.v);
            S = obj.S(obj.v,obj.t);

            m = size(D,1);
            I = speye(m);

            %% Calculate for which indices we need to solve system of equations
            [rows,cols] = find(E);
            j = union(rows,cols);
            i = setdiff(1:m,j);


            %% Calculate matrices need for the timestep
            % Before optimization:  A =  1/k^2 * I - 1/(2*k)*E;
            k = obj.k;
            Aj = 1/k^2 * I(j,j) - 1/(2*k)*E(j,j);
            B =  2/k^2 * I + D;
            C = -1/k^2 * I - 1/(2*k)*E;

            %% Take the timestep
            v = obj.v;

            % Before optimization:  obj.v = A\(B*v + C*v_prev + S);
            obj.v(i) = k^2*(B(i,i)*v(i)   -1/k^2*obj.v_prev(i) + S(i));
            obj.v(j) =  Aj\(B(j,j)*v(j) + C(j,j)*obj.v_prev(j) + S(j));
            obj.v_prev = v;

            %% Update state of the timestepper
            obj.t = obj.t + obj.k;
            obj.n = obj.n + 1;
        end
    end
end