classdef Euler1d < scheme.Scheme
    properties
        m % Number of points in each direction, possibly a vector
        N % Number of points total
        h % Grid spacing
        u % Grid values
        x % Values of x and y for each
        order % Order accuracy for the approximation

        D % non-stabalized scheme operator
        M % Derivative norm
        gamma

        H % Discrete norm
        Hi
        e_l, e_r, e_L, e_R;

    end

    properties (Constant)
        SUBSONIC_INFLOW = 1;
        SUBSONIC_OUTFLOW = -1;
        NO_FLOW = 0;
        SUPERSONIC_INFLOW = 2;
        SUPERSONIC_OUTFLOW = -2;
    end

    methods
        function obj = Euler1d(m,xlim,order,gama,opsGen,do_upwind)
            default_arg('opsGen',@sbp.Ordinary);
            default_arg('gama', 1.4);
            default_arg('do_upwind', false);
            gamma = gama;

            [x, h] = util.get_grid(xlim{:},m);

            if do_upwind
                ops = sbp.Upwind(m,h,order);
                Dp = ops.derivatives.Dp;
                Dm = ops.derivatives.Dm;

                D1 = (Dp + Dm)/2;
                Ddisp = (Dp - Dm)/2;
            else
                ops = opsGen(m,h,order);
                D1 = sparse(ops.derivatives.D1);
            end

            H =  sparse(ops.norms.H);
            Hi = sparse(ops.norms.HI);
            e_l = sparse(ops.boundary.e_1);
            e_r = sparse(ops.boundary.e_m);

            I_x = speye(m);
            I_3 = speye(3);

            D1 = kr(D1, I_3);
            if do_upwind
                Ddisp = kr(Ddisp,I_3);
            end

            % Norms
            obj.H = kr(H,I_3);
            obj.Hi = kr(Hi,I_3);

            % Boundary operators
            obj.e_l  = e_l;
            obj.e_r  = e_r;
            obj.e_L  = kr(e_l,I_3);
            obj.e_R  = kr(e_r,I_3);

            obj.m = m;
            obj.h = h;
            obj.order = order;

            % Man har Q_t+F_x=0 i 1D Euler, där
            % q=[rho, rho*u, e]^T
            % F=[rho*u, rho*u^2+p, (e+p)*u] ^T
            % p=(gamma-1)*(e-rho*u^2/2);


            %Solving on form q_t + F_x = 0

            function o = Fx(q)
                Q = reshape(q,3,m);
                o = reshape(obj.F(Q),3*m,1);
                o = D1*o;
            end

            function o = Fx_disp(q)
                Q = reshape(q,3,m);
                f = reshape(obj.F(Q),3*m,1);

                c = obj.c(Q);
                lambda_max = c+abs(Q(2,:)./Q(1,:));
                % lambda_max = max(lambda_max);

                lamb_Q(1,:) = lambda_max.*Q(1,:);
                lamb_Q(2,:) = lambda_max.*Q(2,:);
                lamb_Q(3,:) = lambda_max.*Q(3,:);

                lamb_q = reshape(lamb_Q,3*m, 1);

                o = -D1*f + Ddisp*lamb_q;
            end

            if do_upwind
                obj.D = @Fx_disp;
            else
                obj.D = @(q)-Fx(q);
            end

            obj.u = x;
            obj.x = kr(x,ones(3,1));
            obj.gamma = gamma;
        end

        % Flux function
        function o = F(obj, Q)
            % Flux: f = [q2; q2.^2/q1 + p(q); (q3+p(q))*q2/q1];
            p = obj.p(Q);
            o = [Q(2,:); Q(2,:).^2./Q(1,:) + p; (Q(3,:)+p).*Q(2,:)./Q(1,:)];
        end

        % Equation of state
        function o = p(obj, Q)
            % Pressure p = (gamma-1)*(q3-q2.^2/q1/2)
            o = (obj.gamma-1)*( Q(3,:)-1/2*Q(2,:).^2./Q(1,:) );
        end

        % Speed of sound
        function o = c(obj, Q)
            % Speed of light c = sqrt(obj.gamma*p/rho);
            o = sqrt(obj.gamma*obj.p(Q)./Q(1,:));
        end

        % Eigen value matrix
        function o = Lambda(obj, q)
            u = q(2)/q(1);
            c = obj.c(q);
            L = [u, u+c, u-c];
            o = diag(L);
        end

        % Diagonalization transformation
        function o = T(obj, q)
            % T is the transformation such that A = T*Lambda*inv(T)
            % where Lambda=diag(u, u+c, u-c)
            rho = q(1);
            u = q(2)/q(1);
            e = q(3);
            gamma = obj.gamma;

            c = sqrt(gamma*obj.p(q)/rho);

            sqrt2gamm = sqrt(2*(gamma-1));

            o = [
                 sqrt2gamm*rho      , rho                               , rho                               ;
                 sqrt2gamm*rho*u    , rho*(u+c)                         , rho*(u-c)                         ;
                 sqrt2gamm*rho*u^2/2, e+(gamma-1)*(e-rho*u^2/2)+rho*u*c , e+(gamma-1)*(e-rho*u^2/2)-rho*u*c ;
            ];
            % Devide columns by stuff to get rid of extra comp?
        end

        function fs = flowStateL(obj, q)
            q_l = obj.e_L'*q;
            c = obj.c(q_l);
            v = q_l(2,:)/q_l(1,:);

            if v > c
                fs = scheme.Euler1d.SUPERSONIC_INFLOW;
            elseif v > 0
                fs = scheme.Euler1d.SUBSONIC_INFLOW;
            elseif v > -c
                fs = scheme.Euler1d.SUBSONIC_OUTFLOW;
            else
                fs = scheme.Euler1d.SUPERSONIC_OUTFLOW;
            end
        end

        % returns positiv values for inlfow, negative for outflow.
        %  +-1 for subsonic
        function fs = flowStateR(obj, q)
            q_r = obj.e_R'*q;
            c = obj.c(q_r);
            v = q_r(2,:)/q_r(1,:);

            if v < -c
                fs = scheme.Euler1d.SUPERSONIC_INFLOW;
            elseif v < 0
                fs = scheme.Euler1d.SUBSONIC_INFLOW;
            elseif v < c
                fs = scheme.Euler1d.SUBSONIC_OUTFLOW;
            else
                fs = scheme.Euler1d.SUPERSONIC_OUTFLOW;
            end
        end

        % Enforces the boundary conditions
        %  w+ = R*w- + g(t)
        function closure = boundary_condition(obj,boundary, type, varargin)
            [e_s,e_S,s] = obj.get_boundary_ops(boundary);

            % Boundary condition on form
            %   w_in = R*w_out + g,       where g is data

            switch type
                case 'wall'
                    closure = obj.boundary_condition_wall(boundary);
                case 'inflow'
                    closure = obj.boundary_condition_inflow(boundary,varargin{:});
                case 'outflow'
                    closure = obj.boundary_condition_outflow(boundary,varargin{:});
                case 'inflow_rho'
                    closure = obj.boundary_condition_inflow_rho(boundary,varargin{:});
                case 'outflow_rho'
                    closure = obj.boundary_condition_outflow_rho(boundary,varargin{:});
                case 'char'
                    closure = obj.boundary_condition_char(boundary,varargin{:});
                otherwise
                    error('Unsupported bc type: %s', type);
            end

        end


        % Sets the boundary condition Lq = g, where
        %   L = L(rho, u, e)
        %   p_in are the indecies of the ingoing characteristics.
        %
        % Returns closure(q,g)
        function closure = boundary_condition_L(obj, boundary, L_fun, p_in)
            [e_s,e_S,s] = obj.get_boundary_ops(boundary);

            p_ot = 1:3;
            p_ot(p_in) = [];

            p = [p_in, p_ot]; % Permutation to sort
            pt(p) = 1:length(p); % Inverse permutation

            function o = closure_fun(q,g)
                % Extract solution at the boundary
                q_s = e_S'*q;
                rho = q_s(1);
                u = q_s(2)/rho;
                e = q_s(3);

                c = obj.c(q_s);

                % Calculate transformation matrix
                T = obj.T(q_s);
                Tin = T(:,p_in);
                Tot = T(:,p_ot);

                % Calculate eigen value matrix
                Lambda = obj.Lambda(q_s);

                % Setup the penalty parameter
                tau1 = -2*abs(Lambda(p_in,p_in));
                tau2 = zeros(length(p_ot),length(p_in)); % Penalty only on ingoing char.

                tauHat = [tau1; tau2];
                tau = e_S*sparse(T*tauHat(pt,:));

                L = L_fun(rho,u,e);

                o = 1/2*obj.Hi * tau * inv(L*Tin)*(L*q_s - g);
            end
            closure = @closure_fun;
        end

        % Return closure(q,g)
        function closure = boundary_condition_char(obj,boundary)
            [e_s,e_S,s] = obj.get_boundary_ops(boundary);

            function o = closure_fun(q, w_data)
                q_s = e_S'*q;
                rho = q_s(1);
                u = q_s(2)/rho;
                e = q_s(3);

                c = obj.c(q_s);

                Lambda = [u, u+c, u-c];

                p_in = find(s*Lambda < 0);
                p_ot = 1:3;
                p_ot(p_in) = [];
                p = [p_in p_ot];
                pt(p) = 1:length(p);

                T = obj.T(q_s);

                tau1 = -2*diag(abs(Lambda(p_in)));
                tau2 = zeros(length(p_ot),length(p_in)); % Penalty only on ingoing char.

                tauHat = [tau1; tau2];

                tau = -s*e_S*sparse(T*tauHat(pt,:));

                w_s = inv(T)*q_s;
                w_in = w_s(p_in);

                w_in_data = w_data(p_in);

                o = 1/2*obj.Hi * tau * (w_in - w_in_data);
            end

            closure = @closure_fun;
        end


        % Return closure(q,[v; p])
        function closure = boundary_condition_inflow(obj, boundary)
            [~,~,s] = obj.get_boundary_ops(boundary);

             switch s
                case -1
                    p_in = [1 2];
                case 1
                    p_in = [1 3];
            end

            a = obj.gamma - 1;
            L = @(rho,u,~)[
                0    1/rho 0;  %v
                0 -1/2*u*a a;  %p
            ];

            closure_raw = boundary_condition_L(obj, boundary, L, g, p_in);
            closure = @(q,p,v) closure_raw(q,[v; p]);
        end

        % Return closure(q, p)
        function closure = boundary_condition_outflow(obj, boundary)
            [~,~,s] = obj.get_boundary_ops(boundary);

             switch s
                case -1
                    p_in = 2;
                case 1
                    p_in = 3;
            end

            a = obj.gamma -1;
            L = @(~,u,~)a*[0 -1/2*u 1];

            closure = boundary_condition_L(obj, boundary, L, p_in);
        end

        % Return closure(q,[v; rho])
        function closure = boundary_condition_inflow_rho(obj, boundary)
            [~,~,s] = obj.get_boundary_ops(boundary);

             switch s
                case -1
                    p_in = [1 2];
                case 1
                    p_in = [1 3];
            end

            a = obj.gamma - 1;
            L = @(rho,~,~)[
                0  1/rho 0;
                1      0 0;
            ];

            closure = boundary_condition_L(obj, boundary, L, p_in);
        end

        % Return closure(q,rho)
        function closure = boundary_condition_outflow_rho(obj, boundary)
            [~,~,s] = obj.get_boundary_ops(boundary);

             switch s
                case -1
                    p_in = 2;
                case 1
                    p_in = 3;
            end

            L = @(~,~,~)[1 0 0];

            closure = boundary_condition_L(obj, boundary, L, p_in);
        end

        % Set wall boundary condition v = 0.
        function closure = boundary_condition_wall(obj,boundary)
            [e_s,e_S,s] = obj.get_boundary_ops(boundary);

            % Vill vi sätta penalty på karateristikan som är nära noll också eller vill
            % vi låta den vara fri?


            switch s
                case -1
                    p_in = 2;
                    p_zero = 1;
                    p_ot = 3;
                case 1
                    p_in = 3;
                    p_zero = 1;
                    p_ot = 2;
                otherwise
                    error();
            end

            p = [p_in, p_zero, p_ot]; % Permutation to sort
            pt(p) = 1:length(p); % Inverse permutation

            function o = closure_fun(q)

                q_s = e_S'*q;
                rho = q_s(1);
                u = q_s(2)/rho;
                c = obj.c(q_s);

                T = obj.T(q_s);
                R = -(u-c)/(u+c);
                % l = [u, u+c, u-c];

                % p_in = find(s*l <= 0);
                % p_ot = find(s*l >  0);


                tau1 = -2*c;
                tau2 = [0; 0]; % Penalty only on ingoing char.

                % Lambda_in = diag(abs(l(p_in)));
                % Lambda_ot = diag(abs(l(p_ot)));

                tauHat = [tau1; tau2];
                tau = -s*e_S*sparse(T*tauHat(pt));

                w_s = inv(T)*q_s;
                % w_s = T\q_s;
                % w_s = Tinv * q_s; % Med analytisk matris
                w_in = w_s(p_in);
                w_ot = w_s(p_ot);

                o = 1/2*obj.Hi * tau * (w_in - R*w_ot);
            end

            closure = @closure_fun;
        end

        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary)
            error('NOT DONE')
            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            [e_u,d1_u,d2_u,d3_u,s_u,gamm_u,delt_u, halfnorm_inv] = obj.get_boundary_ops(boundary);
            [e_v,d1_v,d2_v,d3_v,s_v,gamm_v,delt_v] = neighbour_scheme.get_boundary_ops(neighbour_boundary);

            tuning = 2;

            alpha_u = obj.alpha;
            alpha_v = neighbour_scheme.alpha;

            tau1 = ((alpha_u/2)/delt_u + (alpha_v/2)/delt_v)/2*tuning;
            % tau1 = (alpha_u/2 + alpha_v/2)/(2*delt_u)*tuning;
            tau4 = s_u*alpha_u/2;

            sig2 = ((alpha_u/2)/gamm_u + (alpha_v/2)/gamm_v)/2*tuning;
            sig3 = -s_u*alpha_u/2;

            phi2 = s_u*1/2;

            psi1 = -s_u*1/2;

            tau = tau1*e_u  +                     tau4*d3_u;
            sig =           sig2*d1_u + sig3*d2_u          ;
            phi =           phi2*d1_u                      ;
            psi = psi1*e_u                                 ;

            closure =  halfnorm_inv*(tau*e_u' + sig*d1_u' + phi*alpha_u*d2_u' + psi*alpha_u*d3_u');
            penalty = -halfnorm_inv*(tau*e_v' + sig*d1_v' + phi*alpha_v*d2_v' + psi*alpha_v*d3_v');
        end

        % Ruturns the boundary ops and sign for the boundary specified by the string boundary.
        % The right boundary is considered the positive boundary
        function [e,E,s] = get_boundary_ops(obj,boundary)
            switch boundary
                case 'l'
                    e  = obj.e_l;
                    E  = obj.e_L;
                    s = -1;
                case 'r'
                    e  = obj.e_r;
                    E  = obj.e_R;
                    s = 1;
                otherwise
                    error('No such boundary: boundary = %s',boundary);
            end
        end

        function N = size(obj)
            N = prod(obj.m);
        end

    end
end
