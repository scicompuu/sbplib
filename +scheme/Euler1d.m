classdef SchmBeam2d < noname.Scheme
    properties
        m % Number of points in each direction, possibly a vector
        N % Number of points total
        h % Grid spacing
        u % Grid values
        x % Values of x and y for each
        order % Order accuracy for the approximation

        D % non-stabalized scheme operator
        M % Derivative norm
        alpha

        H % Discrete norm
        Hi
        e_l, e_r

    end

    methods
        function obj = SchmBeam2d(m,xlim,order,gamma,opsGen)
            default_arg('opsGen',@sbp.Ordinary);
            default_arg('gamma', 1.4);

            [x, h] = util.get_grid(xlim{:},m_x);

            ops = opsGen(m_x,h_x,order);

            I_x = speye(m);
            I_3 = speye(3);

            D1 = sparse(ops.derivatives.D1);
            H =  sparse(ops.norms.H);
            Hi = sparse(ops.norms.HI);
            e_l = sparse(ops.boundary.e_1);
            e_r = sparse(ops.boundary.e_m);

            D1 = kr(D1, I_3);

            % Norms
            obj.H = kr(H,I_3);

            % Boundary operators
            obj.e_l  = kr(e_l,I_3);
            obj.e_r  = kr(e_r,I_3);

            obj.m = m;
            obj.h = h;
            obj.order = order;


            % Man har Q_t+F_x=0 i 1D Euler, där
            % q=[rho, rho*u, e]^T
            % F=[rho*u, rho*u^2+p, (e+p)*u] ^T
            % p=(gamma-1)*(e-rho/2*u^2);


            %Solving on form q_t + F_x = 0
            function o = F(q)
                o = [q(2); q(2).^2/q(1) + p(q); (q(3)+p(q))*q(2)/q(1)];
            end

            % Equation of state
            function o = p(q)
                o = (gamma-1)*(q(3)-q(2).^2/q(1)/2);
            end


            % R =
            % [sqrt(2*(gamma-1))*rho      , rho                                , rho           ;
            %  sqrt(2*(gamma-1))*rho*u    , rho*(u+c)                          , rho*(u-c)     ;
            %  sqrt(2*(gamma-1))*rho*u^2/2, e+(gamma-1)*(e-rho*u^2/2)+rho*u*c, e+(gamma-1)*(e-rho*u^2/2)-rho*u*c]);
            function o = R(q)
                rho = q(1);
                u = q(2)/q(1);
                e = q(3);

                sqrt2gamm = sqrt(2*(gamma-1));

                o = [
                     sqrt2gamm*rho      , rho                               , rho                               ;
                     sqrt2gamm*rho*u    , rho*(u+c)                         , rho*(u-c)                         ;
                     sqrt2gamm*rho*u^2/2, e+(gamma-1)*(e-rho*u^2/2)+rho*u*c , e+(gamma-1)*(e-rho*u^2/2)-rho*u*c
                    ];
            end

            function o = Fx(q)
                o = zeros(size(q));
                for i = 1:3:3*m
                    o(i:i+2) = F(q(i:i+2));
                end
            end



            % A=R*Lambda*inv(R), där Lambda=diag(u, u+c, u-c)     (c är ljudhastigheten)
            % c^2=gamma*p/rho
            % function o = A(rho,u,e)
            % end


            obj.D = @Fx;
            obj.u = x;
            obj.x = kr(x,ones(3,1));
        end


        % Closure functions return the opertors applied to the own doamin to close the boundary
        % Penalty functions return the opertors to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of boundary condition if there are several.
        %       data                is a function returning the data that should be applied at the boundary.
        %       neighbour_scheme    is an instance of Scheme that should be interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to interface to.
        function [closure, penalty] = boundary_condition(obj,boundary, alpha,data)
            default_arg('alpha',0);
            default_arg('data',0);

            % Boundary condition on form
            %   w_in = w_out + g,       where g is data

            [e,s] = obj.get_boundary_ops(boundary);

            tuning = 1; % ?????????????????????????

            tau = R(q)*lambda(q)*tuning;     % SHOULD THIS BE abs(lambda)?????

            function closure_fun(q,t)
                q_b = e * q;
            end

            function penalty_fun(q,t)
            end





            % tau1 < -alpha^2/gamma

            tau1 = tuning * alpha/delt;
            tau4 = s*alpha;

            sig2 = tuning * alpha/gamm;
            sig3 = -s*alpha;

            tau = tau1*e+tau4*d3;
            sig = sig2*d1+sig3*d2;

            closure = halfnorm_inv*(tau*e' + sig*d1');

            pp_e = halfnorm_inv*tau;
            pp_d = halfnorm_inv*sig;
            switch class(data)
                case 'double'
                    penalty_e = pp_e*data;
                    penalty_d = pp_d*data;
                case 'function_handle'
                    penalty_e = @(t)pp_e*data(t);
                    penalty_d = @(t)pp_d*data(t);
                otherwise
                    error('Wierd data argument!')
            end

        end

        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary)
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
        function [e,d1,d2,d3,s,gamm, delt, halfnorm_inv] = get_boundary_ops(obj,boundary)
            switch boundary
                case 'w'
                    e  = obj.e_w;
                    d1 = obj.d1_w;
                    d2 = obj.d2_w;
                    d3 = obj.d3_w;
                    s = -1;
                    gamm = obj.gamm_x;
                    delt = obj.delt_x;
                    halfnorm_inv = obj.Hix;
                case 'e'
                    e  = obj.e_e;
                    d1 = obj.d1_e;
                    d2 = obj.d2_e;
                    d3 = obj.d3_e;
                    s = 1;
                    gamm = obj.gamm_x;
                    delt = obj.delt_x;
                    halfnorm_inv = obj.Hix;
                case 's'
                    e  = obj.e_s;
                    d1 = obj.d1_s;
                    d2 = obj.d2_s;
                    d3 = obj.d3_s;
                    s = -1;
                    gamm = obj.gamm_y;
                    delt = obj.delt_y;
                    halfnorm_inv = obj.Hiy;
                case 'n'
                    e  = obj.e_n;
                    d1 = obj.d1_n;
                    d2 = obj.d2_n;
                    d3 = obj.d3_n;
                    s = 1;
                    gamm = obj.gamm_y;
                    delt = obj.delt_y;
                    halfnorm_inv = obj.Hiy;
                otherwise
                    error('No such boundary: boundary = %s',boundary);
            end
        end

        function N = size(obj)
            N = prod(obj.m);
        end

    end
end
