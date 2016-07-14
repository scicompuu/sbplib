classdef Beam < scheme.Scheme
    properties
        order % Order accuracy for the approximation
        grid

        D % non-stabalized scheme operator
        alpha

        H % Discrete norm
        Hi

        e_l, e_r
        d1_l, d1_r
        d2_l, d2_r
        d3_l, d3_r
        gamm
        delt
    end

    methods
        function obj = Beam(grid, order, alpha, opsGen)
            default_arg('alpha', 1);
            % default_arg('opsGen', @sbp.Higher);
            default_arg('opsGen', @sbp.HigherCompatibleVariable); % Supposed to be better

            if ~isa(grid, 'grid.Cartesian') || grid.D() ~= 1
                error('Grid must be 1d cartesian');
            end

            obj.grid = grid;
            obj.order = order;
            obj.alpha = alpha;

            m = grid.m;
            h = grid.scaling();

            ops = opsGen(m, h, order);

            I = speye(m);

            D4 = sparse(ops.derivatives.D4);
            obj.H =  sparse(ops.norms.H);
            obj.Hi = sparse(ops.norms.HI);
            obj.e_l = sparse(ops.boundary.e_1);
            obj.e_r = sparse(ops.boundary.e_m);
            obj.d1_l = sparse(ops.boundary.S_1);
            obj.d1_r = sparse(ops.boundary.S_m);
            obj.d2_l  = sparse(ops.boundary.S2_1);
            obj.d2_r  = sparse(ops.boundary.S2_m);
            obj.d3_l  = sparse(ops.boundary.S3_1);
            obj.d3_r  = sparse(ops.boundary.S3_m);

            obj.D = alpha*D4;

            obj.gamm = h*ops.borrowing.N.S2/2;
            obj.delt = h^3*ops.borrowing.N.S3/2;
        end


        % Closure functions return the opertors applied to the own doamin to close the boundary
        % Penalty functions return the opertors to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of boundary condition if there are several.
        %       neighbour_scheme    is an instance of Scheme that should be interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to interface to.
        function [closure, penalty] = boundary_condition(obj,boundary,type)
            default_arg('type','dn');

            [e, d1, d2, d3, s] = obj.get_boundary_ops(boundary);
            gamm = obj.gamm;
            delt = obj.delt;

            switch type
                case {'dn'} % Dirichlet-neumann boundary condition
                    alpha = obj.alpha;

                    % tau1 < -alpha^2/gamma
                    % tuning = 2;
                    tuning = 1.1;

                    tau1 = tuning * alpha/delt;
                    tau4 = s*alpha;

                    sig2 = tuning * alpha/gamm;
                    sig3 = -s*alpha;

                    tau = tau1*e+tau4*d3;
                    sig = sig2*d1+sig3*d2;

                    closure = obj.Hi*(tau*e' + sig*d1');

                    penalty{1} = -obj.Hi*tau;
                    penalty{2} = -obj.Hi*sig;
                otherwise % Unknown, boundary condition
                    error('No such boundary condition: type = %s',type);
            end
        end

        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary)
            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            [e_u,d1_u,d2_u,d3_u,s_u] = obj.get_boundary_ops(boundary);
            [e_v,d1_v,d2_v,d3_v,s_v] = neighbour_scheme.get_boundary_ops(neighbour_boundary);

            gamm_u = obj.gamm;
            delt_u = obj.delt;

            gamm_v = neighbour_scheme.gamm;
            delt_v = neighbour_scheme.delt;

            % tuning = 2;
            tuning = 1.1;
            % tuning = 0.5;
            % tuning = 0.49998;
            % tuning = 0.3;

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

            closure =  obj.Hi*(tau*e_u' + sig*d1_u' + phi*alpha_u*d2_u' + psi*alpha_u*d3_u');
            penalty = -obj.Hi*(tau*e_v' + sig*d1_v' + phi*alpha_v*d2_v' + psi*alpha_v*d3_v');
        end

        % Returns the boundary ops and sign for the boundary specified by the string boundary.
        % The right boundary is considered the positive boundary
        function [e, d1, d2, d3, s] = get_boundary_ops(obj,boundary)
            switch boundary
                case 'l'
                    e  = obj.e_l;
                    d1 = obj.d1_l;
                    d2 = obj.d2_l;
                    d3 = obj.d3_l;
                    s = -1;
                case 'r'
                    e  = obj.e_r;
                    d1 = obj.d1_r;
                    d2 = obj.d2_r;
                    d3 = obj.d3_r;
                    s = 1;
                otherwise
                    error('No such boundary: boundary = %s',boundary);
            end
        end

        function N = size(obj)
            N = obj.grid.N;
        end

    end
end
