classdef Beam2d < scheme.Scheme
    properties
        grid
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
        d2_w, d2_e, d2_s, d2_n
        d3_w, d3_e, d3_s, d3_n
        gamm_x, gamm_y
        delt_x, delt_y
    end

    methods
        function obj = Beam2d(m,lim,order,alpha,opsGen)
            default_arg('alpha',1);
            default_arg('opsGen',@sbp.Higher);

            if ~isa(grid, 'grid.Cartesian') || grid.D() ~= 2
                error('Grid must be 2d cartesian');
            end

            obj.grid = grid;
            obj.alpha = alpha;
            obj.order = order;

            m_x = grid.m(1);
            m_y = grid.m(2);

            h = grid.scaling();
            h_x = h(1);
            h_y = h(2);

            ops_x = opsGen(m_x,h_x,order);
            ops_y = opsGen(m_y,h_y,order);

            I_x = speye(m_x);
            I_y = speye(m_y);

            D4_x = sparse(ops_x.derivatives.D4);
            H_x =  sparse(ops_x.norms.H);
            Hi_x = sparse(ops_x.norms.HI);
            e_l_x = sparse(ops_x.boundary.e_1);
            e_r_x = sparse(ops_x.boundary.e_m);
            d1_l_x = sparse(ops_x.boundary.S_1);
            d1_r_x = sparse(ops_x.boundary.S_m);
            d2_l_x  = sparse(ops_x.boundary.S2_1);
            d2_r_x  = sparse(ops_x.boundary.S2_m);
            d3_l_x  = sparse(ops_x.boundary.S3_1);
            d3_r_x  = sparse(ops_x.boundary.S3_m);

            D4_y = sparse(ops_y.derivatives.D4);
            H_y =  sparse(ops_y.norms.H);
            Hi_y = sparse(ops_y.norms.HI);
            e_l_y = sparse(ops_y.boundary.e_1);
            e_r_y = sparse(ops_y.boundary.e_m);
            d1_l_y = sparse(ops_y.boundary.S_1);
            d1_r_y = sparse(ops_y.boundary.S_m);
            d2_l_y  = sparse(ops_y.boundary.S2_1);
            d2_r_y  = sparse(ops_y.boundary.S2_m);
            d3_l_y  = sparse(ops_y.boundary.S3_1);
            d3_r_y  = sparse(ops_y.boundary.S3_m);


            D4 = kr(D4_x, I_y) + kr(I_x, D4_y);

            % Norms
            obj.H = kr(H_x,H_y);
            obj.Hx  = kr(H_x,I_x);
            obj.Hy  = kr(I_x,H_y);
            obj.Hix = kr(Hi_x,I_y);
            obj.Hiy = kr(I_x,Hi_y);
            obj.Hi = kr(Hi_x,Hi_y);

            % Boundary operators
            obj.e_w  = kr(e_l_x,I_y);
            obj.e_e  = kr(e_r_x,I_y);
            obj.e_s  = kr(I_x,e_l_y);
            obj.e_n  = kr(I_x,e_r_y);
            obj.d1_w = kr(d1_l_x,I_y);
            obj.d1_e = kr(d1_r_x,I_y);
            obj.d1_s = kr(I_x,d1_l_y);
            obj.d1_n = kr(I_x,d1_r_y);
            obj.d2_w = kr(d2_l_x,I_y);
            obj.d2_e = kr(d2_r_x,I_y);
            obj.d2_s = kr(I_x,d2_l_y);
            obj.d2_n = kr(I_x,d2_r_y);
            obj.d3_w = kr(d3_l_x,I_y);
            obj.d3_e = kr(d3_r_x,I_y);
            obj.d3_s = kr(I_x,d3_l_y);
            obj.d3_n = kr(I_x,d3_r_y);

            obj.D = alpha*D4;

            obj.gamm_x = h_x*ops_x.borrowing.N.S2/2;
            obj.delt_x = h_x^3*ops_x.borrowing.N.S3/2;

            obj.gamm_y = h_y*ops_y.borrowing.N.S2/2;
            obj.delt_y = h_y^3*ops_y.borrowing.N.S3/2;
        end


        % Closure functions return the opertors applied to the own doamin to close the boundary
        % Penalty functions return the opertors to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of boundary condition if there are several.
        %       data                is a function returning the data that should be applied at the boundary.
        %       neighbour_scheme    is an instance of Scheme that should be interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to interface to.
        function [closure, penalty_e,penalty_d] = boundary_condition(obj,boundary,type,data)
            default_arg('type','dn');
            default_arg('data',0);

            [e,d1,d2,d3,s,gamm,delt,halfnorm_inv] = obj.get_boundary_ops(boundary);

            switch type
                % Dirichlet-neumann boundary condition
                case {'dn'}
                    alpha = obj.alpha;

                    % tau1 < -alpha^2/gamma
                    tuning = 1.1;

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

                % Unknown, boundary condition
                otherwise
                    error('No such boundary condition: type = %s',type);
            end
        end

        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary, type)
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
