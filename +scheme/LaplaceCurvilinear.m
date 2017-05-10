classdef LaplaceCurvilinear < scheme.Scheme
    properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing

        grid

        order % Order accuracy for the approximation

        D % non-stabalized scheme operator
        M % Derivative norm
        a,b
        J, Ji
        a11, a12, a22

        H % Discrete norm
        Hi
        H_u, H_v % Norms in the x and y directions
        Hu,Hv % Kroneckerd norms. 1'*Hx*v corresponds to integration in the x dir.
        Hi_u, Hi_v
        Hiu, Hiv
        e_w, e_e, e_s, e_n
        du_w, dv_w
        du_e, dv_e
        du_s, dv_s
        du_n, dv_n
        gamm_u, gamm_v
        lambda

        Dx, Dy % Physical derivatives

        x_u
        x_v
        y_u
        y_v
    end

    methods
        % Implements  a*div(b*grad(u)) as a SBP scheme
        function obj = LaplaceCurvilinear(g ,order, a, b, opSet)
            default_arg('opSet',@sbp.D2Variable);
            default_arg('a', 1);
            default_arg('b', 1);


            if b ~=1
                error('Not implemented yet')
            end

            assert(isa(g, 'grid.Curvilinear'))

            m = g.size();
            m_u = m(1);
            m_v = m(2);
            m_tot = g.N();

            h = g.scaling();
            h_u = h(1);
            h_v = h(2);

            % Operators
            ops_u = opSet(m_u, {0, 1}, order);
            ops_v = opSet(m_v, {0, 1}, order);

            I_u = speye(m_u);
            I_v = speye(m_v);

            D1_u = ops_u.D1;
            D2_u = ops_u.D2;
            H_u =  ops_u.H;
            Hi_u = ops_u.HI;
            e_l_u = ops_u.e_l;
            e_r_u = ops_u.e_r;
            d1_l_u = ops_u.d1_l;
            d1_r_u = ops_u.d1_r;

            D1_v = ops_v.D1;
            D2_v = ops_v.D2;
            H_v =  ops_v.H;
            Hi_v = ops_v.HI;
            e_l_v = ops_v.e_l;
            e_r_v = ops_v.e_r;
            d1_l_v = ops_v.d1_l;
            d1_r_v = ops_v.d1_r;

            Du = kr(D1_u,I_v);
            Dv = kr(I_u,D1_v);

            % Metric derivatives
            coords = g.points();
            x = coords(:,1);
            y = coords(:,2);

            x_u = Du*x;
            x_v = Dv*x;
            y_u = Du*y;
            y_v = Dv*y;

            J = x_u.*y_v - x_v.*y_u;
            a11 =  1./J .* (x_v.^2  + y_v.^2);
            a12 = -1./J .* (x_u.*x_v + y_u.*y_v);
            a22 =  1./J .* (x_u.^2  + y_u.^2);
            lambda = 1/2 * (a11 + a22 - sqrt((a11-a22).^2 + 4*a12.^2));

            % Assemble full operators
            L_12 = spdiags(a12, 0, m_tot, m_tot);
            Duv = Du*L_12*Dv;
            Dvu = Dv*L_12*Du;

            Duu = sparse(m_tot);
            Dvv = sparse(m_tot);
            ind = grid.funcToMatrix(g, 1:m_tot);

            for i = 1:m_v
                D = D2_u(a11(ind(:,i)));
                p = ind(:,i);
                Duu(p,p) = D;
            end

            for i = 1:m_u
                D = D2_v(a22(ind(i,:)));
                p = ind(i,:);
                Dvv(p,p) = D;
            end

            obj.H = kr(H_u,H_v);
            obj.Hi = kr(Hi_u,Hi_v);
            obj.Hu  = kr(H_u,I_v);
            obj.Hv  = kr(I_u,H_v);
            obj.Hiu = kr(Hi_u,I_v);
            obj.Hiv = kr(I_u,Hi_v);

            obj.e_w  = kr(e_l_u,I_v);
            obj.e_e  = kr(e_r_u,I_v);
            obj.e_s  = kr(I_u,e_l_v);
            obj.e_n  = kr(I_u,e_r_v);
            obj.du_w = kr(d1_l_u,I_v);
            obj.dv_w = (obj.e_w'*Dv)';
            obj.du_e = kr(d1_r_u,I_v);
            obj.dv_e = (obj.e_e'*Dv)';
            obj.du_s = (obj.e_s'*Du)';
            obj.dv_s = kr(I_u,d1_l_v);
            obj.du_n = (obj.e_n'*Du)';
            obj.dv_n = kr(I_u,d1_r_v);

            obj.x_u = x_u;
            obj.x_v = x_v;
            obj.y_u = y_u;
            obj.y_v = y_v;

            obj.m = m;
            obj.h = [h_u h_v];
            obj.order = order;
            obj.grid = g;

            obj.a = a;
            obj.b = b;
            obj.J = spdiags(J, 0, m_tot, m_tot);
            obj.Ji = spdiags(1./J, 0, m_tot, m_tot);
            obj.a11 = a11;
            obj.a12 = a12;
            obj.a22 = a22;
            obj.D = obj.Ji*a*(Duu + Duv + Dvu + Dvv);
            obj.lambda = lambda;

            obj.Dx = spdiag( y_v./J)*Du + spdiag(-y_u./J)*Dv;
            obj.Dy = spdiag(-x_v./J)*Du + spdiag( x_u./J)*Dv;

            obj.gamm_u = h_u*ops_u.borrowing.M.d1;
            obj.gamm_v = h_v*ops_v.borrowing.M.d1;
        end


        % Closure functions return the opertors applied to the own doamin to close the boundary
        % Penalty functions return the opertors to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of boundary condition if there are several.
        %       data                is a function returning the data that should be applied at the boundary.
        %       neighbour_scheme    is an instance of Scheme that should be interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to interface to.
        function [closure, penalty] = boundary_condition(obj, boundary, type, parameter)
            default_arg('type','neumann');
            default_arg('parameter', []);

            [e, d_n, d_t, coeff_n, coeff_t, s, gamm, halfnorm_inv  ,              ~,          ~, ~, scale_factor] = obj.get_boundary_ops(boundary);
            switch type
                % Dirichlet boundary condition
                case {'D','d','dirichlet'}
                    % v denotes the solution in the neighbour domain
                    tuning = 1.2;
                    % tuning = 20.2;
                    [e, d_n, d_t, coeff_n, coeff_t, s, gamm, halfnorm_inv_n, halfnorm_inv_t, halfnorm_t] = obj.get_boundary_ops(boundary);

                    a_n = spdiag(coeff_n);
                    a_t = spdiag(coeff_t);

                    F = (s * a_n * d_n' + s * a_t*d_t')';

                    u = obj;

                    b1 = gamm*u.lambda./u.a11.^2;
                    b2 = gamm*u.lambda./u.a22.^2;

                    tau  = -1./b1 - 1./b2;
                    tau = tuning * spdiag(tau);
                    sig1 = 1;

                    penalty_parameter_1 = halfnorm_inv_n*(tau + sig1*halfnorm_inv_t*F*e'*halfnorm_t)*e;

                    closure = obj.Ji*obj.a * penalty_parameter_1*e';
                    penalty = -obj.Ji*obj.a * penalty_parameter_1;


                % Neumann boundary condition
                case {'N','n','neumann'}
                    a_n = spdiags(coeff_n,0,length(coeff_n),length(coeff_n));
                    a_t = spdiags(coeff_t,0,length(coeff_t),length(coeff_t));
                    d = (a_n * d_n' + a_t*d_t')';

                    tau1 = -s;
                    tau2 = 0;
                    tau = obj.a * obj.Ji*(tau1*e + tau2*d);

                    closure = halfnorm_inv*tau*d';
                    penalty = -halfnorm_inv*tau;

                % Characteristic boundary condition
                case {'characteristic', 'char', 'c'}
                    default_arg('parameter', 1);
                    beta = parameter;

                    a_n = spdiags(coeff_n,0,length(coeff_n),length(coeff_n));
                    a_t = spdiags(coeff_t,0,length(coeff_t),length(coeff_t));
                    d = s*(a_n * d_n' + a_t*d_t')'; % outward facing normal derivative

                    tau = -obj.a * 1/beta*obj.Ji*e;

                    closure{1} = halfnorm_inv*tau*spdiag(scale_factor)*e';
                    closure{2} = halfnorm_inv*tau*beta*d';
                    penalty = -halfnorm_inv*tau;

                % Unknown, boundary condition
                otherwise
                    error('No such boundary condition: type = %s',type);
            end
        end

        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary)
            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            tuning = 1.2;
            % tuning = 20.2;
            [e_u, d_n_u, d_t_u, coeff_n_u, coeff_t_u, s_u, gamm_u, halfnorm_inv_u_n, halfnorm_inv_u_t, halfnorm_u_t, I_u] = obj.get_boundary_ops(boundary);
            [e_v, d_n_v, d_t_v, coeff_n_v, coeff_t_v, s_v, gamm_v, halfnorm_inv_v_n, halfnorm_inv_v_t, halfnorm_v_t, I_v] = neighbour_scheme.get_boundary_ops(neighbour_boundary);

            a_n_u = spdiag(coeff_n_u);
            a_t_u = spdiag(coeff_t_u);
            a_n_v = spdiag(coeff_n_v);
            a_t_v = spdiag(coeff_t_v);

            F_u = (s_u * a_n_u * d_n_u' + s_u * a_t_u*d_t_u')';
            F_v = (s_v * a_n_v * d_n_v' + s_v * a_t_v*d_t_v')';

            u = obj;
            v = neighbour_scheme;

            b1_u = gamm_u*u.lambda(I_u)./u.a11(I_u).^2;
            b2_u = gamm_u*u.lambda(I_u)./u.a22(I_u).^2;
            b1_v = gamm_v*v.lambda(I_v)./v.a11(I_v).^2;
            b2_v = gamm_v*v.lambda(I_v)./v.a22(I_v).^2;

            tau = -1./(4*b1_u) -1./(4*b1_v) -1./(4*b2_u) -1./(4*b2_v);
            tau = tuning * spdiag(tau);
            sig1 = 1/2;
            sig2 = -1/2;

            penalty_parameter_1 = halfnorm_inv_u_n*(e_u*tau + sig1*halfnorm_inv_u_t*F_u*e_u'*halfnorm_u_t*e_u);
            penalty_parameter_2 = halfnorm_inv_u_n * sig2 * e_u;


            closure = obj.Ji*obj.a * ( penalty_parameter_1*e_u' + penalty_parameter_2*F_u');
            penalty = obj.Ji*obj.a * (-penalty_parameter_1*e_v' + penalty_parameter_2*F_v');
        end

        % Ruturns the boundary ops and sign for the boundary specified by the string boundary.
        % The right boundary is considered the positive boundary
        %
        %  I -- the indecies of the boundary points in the grid matrix
        function [e, d_n, d_t, coeff_n, coeff_t, s, gamm, halfnorm_inv_n, halfnorm_inv_t, halfnorm_t, I, scale_factor] = get_boundary_ops(obj, boundary)

            % gridMatrix = zeros(obj.m(2),obj.m(1));
            % gridMatrix(:) = 1:numel(gridMatrix);

            ind = grid.funcToMatrix(obj.grid, 1:prod(obj.m));

            switch boundary
                case 'w'
                    e = obj.e_w;
                    d_n = obj.du_w;
                    d_t = obj.dv_w;
                    s = -1;

                    I = ind(1,:);
                    coeff_n = obj.a11(I);
                    coeff_t = obj.a12(I);
                    scale_factor = sqrt(obj.x_v(I).^2 + obj.y_v(I).^2);
                case 'e'
                    e = obj.e_e;
                    d_n = obj.du_e;
                    d_t = obj.dv_e;
                    s = 1;

                    I = ind(end,:);
                    coeff_n = obj.a11(I);
                    coeff_t = obj.a12(I);
                    scale_factor = sqrt(obj.x_v(I).^2 + obj.y_v(I).^2);
                case 's'
                    e = obj.e_s;
                    d_n = obj.dv_s;
                    d_t = obj.du_s;
                    s = -1;

                    I = ind(:,1)';
                    coeff_n = obj.a22(I);
                    coeff_t = obj.a12(I);
                    scale_factor = sqrt(obj.x_u(I).^2 + obj.y_u(I).^2);
                case 'n'
                    e = obj.e_n;
                    d_n = obj.dv_n;
                    d_t = obj.du_n;
                    s = 1;

                    I = ind(:,end)';
                    coeff_n = obj.a22(I);
                    coeff_t = obj.a12(I);
                    scale_factor = sqrt(obj.x_u(I).^2 + obj.y_u(I).^2);
                otherwise
                    error('No such boundary: boundary = %s',boundary);
            end

            switch boundary
                case {'w','e'}
                    halfnorm_inv_n = obj.Hiu;
                    halfnorm_inv_t = obj.Hiv;
                    halfnorm_t = obj.Hv;
                    gamm = obj.gamm_u;
                case {'s','n'}
                    halfnorm_inv_n = obj.Hiv;
                    halfnorm_inv_t = obj.Hiu;
                    halfnorm_t = obj.Hu;
                    gamm = obj.gamm_v;
            end
        end

        function N = size(obj)
            N = prod(obj.m);
        end


    end
end