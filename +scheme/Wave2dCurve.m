classdef Wave2dCurve < scheme.Scheme
    properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing
        u,v % Grid
        x,y % Values of x and y for each grid point
        X,Y % Grid point locations as matrices
        order % Order accuracy for the approximation

        D % non-stabalized scheme operator
        M % Derivative norm
        c
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
    end

    methods
        function obj = Wave2dCurve(m,ti,order,c,opSet)
            default_arg('opSet',@sbp.Variable);
            default_arg('c', 1);

            if length(m) == 1
                m = [m m];
            end

            m_u = m(1);
            m_v = m(2);
            m_tot = m_u*m_v;

            [u, h_u] = util.get_grid(0, 1, m_u);
            [v, h_v] = util.get_grid(0, 1, m_v);


            % Operators
            ops_u = opSet(m_u,h_u,order);
            ops_v = opSet(m_v,h_v,order);

            I_u = speye(m_u);
            I_v = speye(m_v);

            D1_u = sparse(ops_u.derivatives.D1);
            D2_u = ops_u.derivatives.D2;
            H_u =  sparse(ops_u.norms.H);
            Hi_u = sparse(ops_u.norms.HI);
            % M_u =  sparse(ops_u.norms.M);
            e_l_u = sparse(ops_u.boundary.e_1);
            e_r_u = sparse(ops_u.boundary.e_m);
            d1_l_u = sparse(ops_u.boundary.S_1);
            d1_r_u = sparse(ops_u.boundary.S_m);

            D1_v = sparse(ops_v.derivatives.D1);
            D2_v = ops_v.derivatives.D2;
            H_v =  sparse(ops_v.norms.H);
            Hi_v = sparse(ops_v.norms.HI);
            % M_v =  sparse(ops_v.norms.M);
            e_l_v = sparse(ops_v.boundary.e_1);
            e_r_v = sparse(ops_v.boundary.e_m);
            d1_l_v = sparse(ops_v.boundary.S_1);
            d1_r_v = sparse(ops_v.boundary.S_m);


            % Metric derivatives
            [X,Y] = ti.map(u,v);

            [x_u,x_v] = gridDerivatives(X,D1_u,D1_v);
            [y_u,y_v] = gridDerivatives(Y,D1_u,D1_v);



            J = x_u.*y_v - x_v.*y_u;
            a11 =  1./J .* (x_v.^2  + y_v.^2);  %% GÖR SOM MATRISER
            a12 = -1./J .* (x_u.*x_v + y_u.*y_v);
            a22 =  1./J .* (x_u.^2  + y_u.^2);
            lambda = 1/2 * (a11 + a22 - sqrt((a11-a22).^2 + 4*a12.^2));

            dof_order = reshape(1:m_u*m_v,m_v,m_u);

            Duu = sparse(m_tot);
            Dvv = sparse(m_tot);

            for i = 1:m_v
                D = D2_u(a11(i,:));
                p = dof_order(i,:);
                Duu(p,p) = D;
            end

            for i = 1:m_u
                D = D2_v(a22(:,i));
                p = dof_order(:,i);
                Dvv(p,p) = D;
            end

            L_12 = spdiags(a12(:),0,m_tot,m_tot);
            Du = kr(D1_u,I_v);
            Dv = kr(I_u,D1_v);

            Duv = Du*L_12*Dv;
            Dvu = Dv*L_12*Du;



            obj.H = kr(H_u,H_v);
            obj.Hi = kr(Hi_u,Hi_v);
            obj.Hu  = kr(H_u,I_v);
            obj.Hv  = kr(I_u,H_v);
            obj.Hiu = kr(Hi_u,I_v);
            obj.Hiv = kr(I_u,Hi_v);

            % obj.M = kr(M_u,H_v)+kr(H_u,M_v);
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

            obj.m = m;
            obj.h = [h_u h_v];
            obj.order = order;


            obj.c = c;
            obj.J = spdiags(J(:),0,m_tot,m_tot);
            obj.Ji = spdiags(1./J(:),0,m_tot,m_tot);
            obj.a11 = a11;
            obj.a12 = a12;
            obj.a22 = a22;
            obj.D = obj.Ji*c^2*(Duu + Duv + Dvu + Dvv);
            obj.u = u;
            obj.v = v;
            obj.X = X;
            obj.Y = Y;
            obj.x = X(:);
            obj.y = Y(:);
            obj.lambda = lambda;

            obj.gamm_u = h_u*ops_u.borrowing.M.S;
            obj.gamm_v = h_v*ops_v.borrowing.M.S;
        end


        % Closure functions return the opertors applied to the own doamin to close the boundary
        % Penalty functions return the opertors to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of boundary condition if there are several.
        %       data                is a function returning the data that should be applied at the boundary.
        %       neighbour_scheme    is an instance of Scheme that should be interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to interface to.
        function [closure, penalty] = boundary_condition(obj,boundary,type,data)
            default_arg('type','neumann');
            default_arg('data',0);

            [e, d_n, d_t, coeff_n, coeff_t, s, gamm, halfnorm_inv] = obj.get_boundary_ops(boundary);

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
                    tau = tuning * spdiag(tau(:));
                    sig1 = 1/2;

                    penalty_parameter_1 = halfnorm_inv_n*(tau + sig1*halfnorm_inv_t*F*e'*halfnorm_t)*e;

                    closure = obj.Ji*obj.c^2 * penalty_parameter_1*e';
                    pp = -obj.Ji*obj.c^2 * penalty_parameter_1;
                    switch class(data)
                        case 'double'
                            penalty = pp*data;
                        case 'function_handle'
                            penalty = @(t)pp*data(t);
                        otherwise
                            error('Weird data argument!')
                    end


                % Neumann boundary condition
                case {'N','n','neumann'}
                    c = obj.c;


                    a_n = spdiags(coeff_n,0,length(coeff_n),length(coeff_n));
                    a_t = spdiags(coeff_t,0,length(coeff_t),length(coeff_t));
                    d = (a_n * d_n' + a_t*d_t')';

                    tau1 = -s;
                    tau2 = 0;
                    tau = c.^2 * obj.Ji*(tau1*e + tau2*d);

                    closure = halfnorm_inv*tau*d';

                    pp = halfnorm_inv*tau;
                    switch class(data)
                        case 'double'
                            penalty = pp*data;
                        case 'function_handle'
                            penalty = @(t)pp*data(t);
                        otherwise
                            error('Weird data argument!')
                    end

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
            tau = tuning * spdiag(tau(:)); % Probably correct until here, see eq 27
            sig1 = 1/2;
            sig2 = -1/2;

            % penalty_parameter_1 = halfnorm_inv_u_n*(tau + sig1*halfnorm_inv_u_t*F_u*e_u'*halfnorm_u_t)*e_u;  %% This is what is in the paper, but there is an error in dimensions.
            penalty_parameter_1 = halfnorm_inv_u_n*(e_u*tau + sig1*halfnorm_inv_u_t*F_u*e_u'*halfnorm_u_t*e_u); %% Random guess at a fix, should check theory for this.
            penalty_parameter_2 = halfnorm_inv_u_n * sig2 * e_u;


            closure = obj.Ji*obj.c^2 * ( penalty_parameter_1*e_u' + penalty_parameter_2*F_u');
            penalty = obj.Ji*obj.c^2 * (-penalty_parameter_1*e_v' + penalty_parameter_2*F_v');
        end

        % Ruturns the boundary ops and sign for the boundary specified by the string boundary.
        % The right boundary is considered the positive boundary
        %
        %  I -- the indecies of the boundary points in the grid matrix
        function [e, d_n, d_t, coeff_n, coeff_t, s, gamm, halfnorm_inv_n, halfnorm_inv_t, halfnorm_t, I] = get_boundary_ops(obj,boundary)

            gridMatrix = zeros(obj.m(2),obj.m(1));
            gridMatrix(:) = 1:numel(gridMatrix);

            switch boundary
                case 'w'
                    e = obj.e_w;
                    d_n = obj.du_w;
                    d_t = obj.dv_w;
                    s = -1;

                    I = gridMatrix(:,1);
                    coeff_n = obj.a11(:,1);
                    coeff_t = obj.a12(:,1);
                case 'e'
                    e = obj.e_e;
                    d_n = obj.du_e;
                    d_t = obj.dv_e;
                    s = 1;

                    I = gridMatrix(:,end);
                    coeff_n = obj.a11(:,end);
                    coeff_t = obj.a12(:,end);
                case 's'
                    e = obj.e_s;
                    d_n = obj.dv_s;
                    d_t = obj.du_s;
                    s = -1;

                    I = gridMatrix(1,:)';
                    coeff_n = obj.a22(1,:)';
                    coeff_t = obj.a12(1,:)';
                case 'n'
                    e = obj.e_n;
                    d_n = obj.dv_n;
                    d_t = obj.du_n;
                    s = 1;

                    I = gridMatrix(end,:)';
                    coeff_n = obj.a22(end,:)';
                    coeff_t = obj.a12(end,:)';
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

    methods(Static)
        % Calculates the matrcis need for the inteface coupling between boundary bound_u of scheme schm_u
        % and bound_v of scheme schm_v.
        %   [uu, uv, vv, vu] = inteface_couplong(A,'r',B,'l')
        function [uu, uv, vv, vu] = interface_coupling(schm_u,bound_u,schm_v,bound_v)
            [uu,uv] = schm_u.interface(bound_u,schm_v,bound_v);
            [vv,vu] = schm_v.interface(bound_v,schm_u,bound_u);
        end
    end
end