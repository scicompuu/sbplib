classdef Beam < scheme.Scheme
    properties
        order % Order accuracy for the approximation
        grid

        D % non-stabalized scheme operator
        alpha

        h
        H % Discrete norm
        Hi

        e_l,  e_r
        d1_l, d1_r
        d2_l, d2_r
        d3_l, d3_r
        gamm
        delt
        alphaII
        alphaIII

        opt
    end

    methods
        function obj = Beam(grid, order, alpha, opsGen, opt)
            default_arg('alpha', -1);

            % default_arg('opsGen', @sbp.D4);
            default_arg('opsGen', @sbp.D4Variable); % Supposed to be better

            opt_default.interface_l.tuning = 1.1;
            opt_default.interface_l.tau = [];
            opt_default.interface_l.sig = [];
            opt_default.interface_r.tuning = 1.1;
            opt_default.interface_r.tau = [];
            opt_default.interface_r.sig = [];
            default_struct('opt', opt_default);

            if ~isa(grid, 'grid.Cartesian') || grid.D() ~= 1
                error('Grid must be 1d cartesian');
            end

            obj.grid = grid;
            obj.order = order;
            obj.alpha = alpha;

            m = grid.m;
            h = grid.scaling();

            x_lim = {grid.x{1}(1), grid.x{1}(end)};
            ops = opsGen(m, x_lim, order);

            D4       = ops.D4;
            obj.H    = ops.H;
            obj.Hi   = ops.HI;
            obj.e_l  = ops.e_l;
            obj.e_r  = ops.e_r;
            obj.d1_l = ops.d1_l;
            obj.d1_r = ops.d1_r;
            obj.d2_l = ops.d2_l;
            obj.d2_r = ops.d2_r;
            obj.d3_l = ops.d3_l;
            obj.d3_r = ops.d3_r;

            obj.D = alpha*D4;

            alphaII  = ops.borrowing.N.S2/2;
            alphaIII = ops.borrowing.N.S3/2;

            obj.gamm = h*alphaII;
            obj.delt = h^3*alphaIII;
            obj.alphaII = alphaII;
            obj.alphaIII = alphaIII;
            obj.h = h;
            obj.opt = opt;
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
                case {'dn', 'clamped'} % Dirichlet-neumann boundary condition
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


                case {'free'}
                    a = obj.alpha;

                    tau =  s*a*d1;
                    sig = -s*a*e;

                    closure = obj.Hi*(tau*d2' + sig*d3');
                    penalty{1} = -obj.Hi*tau;
                    penalty{1} = -obj.Hi*sig;


                otherwise % Unknown, boundary condition
                    error('No such boundary condition: type = %s',type);
            end
        end

        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary)
            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            [e_u,d1_u,d2_u,d3_u,s_u] = obj.get_boundary_ops(boundary);
            [e_v,d1_v,d2_v,d3_v,s_v] = neighbour_scheme.get_boundary_ops(neighbour_boundary);


            alpha_u = obj.alpha;
            alpha_v = neighbour_scheme.alpha;


            switch boundary
                case 'l'
                    interface_opt = obj.opt.interface_l;
                case 'r'
                    interface_opt = obj.opt.interface_r;
            end


            if isempty(interface_opt.tau) && isempty(interface_opt.sig)
                gamm_u = obj.gamm;
                delt_u = obj.delt;

                gamm_v = neighbour_scheme.gamm;
                delt_v = neighbour_scheme.delt;

                tuning = interface_opt.tuning;

                tau1 = ((alpha_u/2)/delt_u + (alpha_v/2)/delt_v)/2*tuning;
                sig2 = ((alpha_u/2)/gamm_u + (alpha_v/2)/gamm_v)/2*tuning;
            else
                h_u = obj.h;
                h_v = neighbour_scheme.h;

                switch neighbour_boundary
                    case 'l'
                        neighbour_interface_opt = neighbour_scheme.opt.interface_l;
                    case 'r'
                        neighbour_interface_opt = neighbour_scheme.opt.interface_r;
                end

                tau_u = interface_opt.tau;
                sig_u = interface_opt.sig;
                tau_v = neighbour_interface_opt.tau;
                sig_v = neighbour_interface_opt.sig;

                tau1 = tau_u/h_u^3 + tau_v/h_v^3;
                sig2 = sig_u/h_u   + sig_v/h_v;
            end

            tau4 = s_u*alpha_u/2;
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
