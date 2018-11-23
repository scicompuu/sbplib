classdef Wave2d < scheme.Scheme
    properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing
        x,y % Grid
        X,Y % Values of x and y for each grid point
        order % Order accuracy for the approximation

        D % non-stabalized scheme operator
        M % Derivative norm
        alpha

        H % Discrete norm
        Hi
        H_x, H_y % Norms in the x and y directions
        Hx,Hy % Kroneckerd norms. 1'*Hx*v corresponds to integration in the x dir.
        Hi_x, Hi_y
        Hix, Hiy
        e_w, e_e, e_s, e_n
        d1_w, d1_e, d1_s, d1_n
        gamm_x, gamm_y
    end

    methods
        function obj = Wave2d(m,lim,order,alpha)
            default_arg('alpha',1);

            xlim = lim{1};
            ylim = lim{2};

            if length(m) == 1
                m = [m m];
            end

            m_x = m(1);
            m_y = m(2);

            [x, h_x] = util.get_grid(xlim{:},m_x);
            [y, h_y] = util.get_grid(ylim{:},m_y);

            ops_x = sbp.Ordinary(m_x,h_x,order);
            ops_y = sbp.Ordinary(m_y,h_y,order);

            I_x = speye(m_x);
            I_y = speye(m_y);

            D2_x = sparse(ops_x.derivatives.D2);
            H_x =  sparse(ops_x.norms.H);
            Hi_x = sparse(ops_x.norms.HI);
            M_x =  sparse(ops_x.norms.M);
            e_l_x = sparse(ops_x.boundary.e_1);
            e_r_x = sparse(ops_x.boundary.e_m);
            d1_l_x = sparse(ops_x.boundary.S_1);
            d1_r_x = sparse(ops_x.boundary.S_m);

            D2_y = sparse(ops_y.derivatives.D2);
            H_y =  sparse(ops_y.norms.H);
            Hi_y = sparse(ops_y.norms.HI);
            M_y =  sparse(ops_y.norms.M);
            e_l_y = sparse(ops_y.boundary.e_1);
            e_r_y = sparse(ops_y.boundary.e_m);
            d1_l_y = sparse(ops_y.boundary.S_1);
            d1_r_y = sparse(ops_y.boundary.S_m);

            D2 = kr(D2_x, I_y) + kr(I_x, D2_y);
            obj.H = kr(H_x,H_y);
            obj.Hx  = kr(H_x,I_y);
            obj.Hy  = kr(I_x,H_y);
            obj.Hix = kr(Hi_x,I_y);
            obj.Hiy = kr(I_x,Hi_y);
            obj.Hi = kr(Hi_x,Hi_y);
            obj.M = kr(M_x,H_y)+kr(H_x,M_y);
            obj.e_w  = kr(e_l_x,I_y);
            obj.e_e  = kr(e_r_x,I_y);
            obj.e_s  = kr(I_x,e_l_y);
            obj.e_n  = kr(I_x,e_r_y);
            obj.d1_w = kr(d1_l_x,I_y);
            obj.d1_e = kr(d1_r_x,I_y);
            obj.d1_s = kr(I_x,d1_l_y);
            obj.d1_n = kr(I_x,d1_r_y);

            obj.m = m;
            obj.h = [h_x h_y];
            obj.order = order;

            obj.alpha = alpha;
            obj.D = alpha*D2;
            obj.x = x;
            obj.y = y;
            obj.X = kr(x,ones(m_y,1));
            obj.Y = kr(ones(m_x,1),y);

            obj.gamm_x = h_x*ops_x.borrowing.M.S;
            obj.gamm_y = h_y*ops_y.borrowing.M.S;
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

            [e,d,s,gamm,halfnorm_inv] = obj.get_boundary_ops(boundary);

            switch type
                % Dirichlet boundary condition
                case {'D','d','dirichlet'}
                    alpha = obj.alpha;

                    % tau1 < -alpha^2/gamma
                    tuning = 1.1;
                    tau1 = -tuning*alpha/gamm;
                    tau2 =  s*alpha;

                    p = tau1*e + tau2*d;

                    closure = halfnorm_inv*p*e';

                    pp = halfnorm_inv*p;
                    switch class(data)
                        case 'double'
                            penalty = pp*data;
                        case 'function_handle'
                            penalty = @(t)pp*data(t);
                        otherwise
                            error('Wierd data argument!')
                    end


                % Neumann boundary condition
                case {'N','n','neumann'}
                    alpha = obj.alpha;
                    tau1 = -s*alpha;
                    tau2 = 0;
                    tau = tau1*e + tau2*d;

                    closure = halfnorm_inv*tau*d';

                    pp = halfnorm_inv*tau;
                    switch class(data)
                        case 'double'
                            penalty = pp*data;
                        case 'function_handle'
                            penalty = @(t)pp*data(t);
                        otherwise
                            error('Wierd data argument!')
                    end

                % Unknown, boundary condition
                otherwise
                    error('No such boundary condition: type = %s',type);
            end
        end

        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary,type)
            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            [e_u,d_u,s_u,gamm_u, halfnorm_inv] = obj.get_boundary_ops(boundary);
            [e_v,d_v,s_v,gamm_v] = neighbour_scheme.get_boundary_ops(neighbour_boundary);

            tuning = 1.1;

            alpha_u = obj.alpha;
            alpha_v = neighbour_scheme.alpha;

            % tau1 < -(alpha_u/gamm_u + alpha_v/gamm_v)

            tau1 = -(alpha_u/gamm_u + alpha_v/gamm_v) * tuning;
            tau2 = s_u*1/2*alpha_u;
            sig1 = s_u*(-1/2);
            sig2 = 0;

            tau = tau1*e_u + tau2*d_u;
            sig = sig1*e_u + sig2*d_u;

            closure = halfnorm_inv*( tau*e_u' + sig*alpha_u*d_u');
            penalty = halfnorm_inv*(-tau*e_v' - sig*alpha_v*d_v');
        end

        % Ruturns the boundary ops and sign for the boundary specified by the string boundary.
        % The right boundary is considered the positive boundary
        function [e,d,s,gamm, halfnorm_inv] = get_boundary_ops(obj,boundary)
            switch boundary
                case 'w'
                    e = obj.e_w;
                    d = obj.d1_w;
                    s = -1;
                    gamm = obj.gamm_x;
                    halfnorm_inv = obj.Hix;
                case 'e'
                    e = obj.e_e;
                    d = obj.d1_e;
                    s = 1;
                    gamm = obj.gamm_x;
                    halfnorm_inv = obj.Hix;
                case 's'
                    e = obj.e_s;
                    d = obj.d1_s;
                    s = -1;
                    gamm = obj.gamm_y;
                    halfnorm_inv = obj.Hiy;
                case 'n'
                    e = obj.e_n;
                    d = obj.d1_n;
                    s = 1;
                    gamm = obj.gamm_y;
                    halfnorm_inv = obj.Hiy;
                otherwise
                    error('No such boundary: boundary = %s',boundary);
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