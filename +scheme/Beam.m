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

        opt % TODO: Get rid of this and use the interface type instead
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

            e  = obj.getBoundaryOperator('e',  boundary);
            d1 = obj.getBoundaryOperator('d1', boundary);
            d2 = obj.getBoundaryOperator('d2', boundary);
            d3 = obj.getBoundaryOperator('d3', boundary);
            s = obj.getBoundarySign(boundary);
            gamm = obj.gamm;
            delt = obj.delt;


            % TODO: Can this be simplifed? Can I handle conditions on u on its own, u_x on its own ...

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
                    penalty{2} = -obj.Hi*sig;

                case 'e'
                    alpha = obj.alpha;
                    tuning = 1.1;

                    tau1 = tuning * alpha/delt;
                    tau4 = s*alpha;

                    tau = tau1*e+tau4*d3;

                    closure = obj.Hi*tau*e';
                    penalty = -obj.Hi*tau;
                case 'd1'
                    alpha = obj.alpha;

                    tuning = 1.1;

                    sig2 = tuning * alpha/gamm;
                    sig3 = -s*alpha;

                    sig = sig2*d1+sig3*d2;

                    closure = obj.Hi*sig*d1';
                    penalty = -obj.Hi*sig;

                case 'd2'
                    a = obj.alpha;

                    tau =  s*a*d1;

                    closure = obj.Hi*tau*d2';
                    penalty = -obj.Hi*tau;
                case 'd3'
                    a = obj.alpha;

                    sig = -s*a*e;

                    closure = obj.Hi*sig*d3';
                    penalty = -obj.Hi*sig;

                otherwise % Unknown, boundary condition
                    error('No such boundary condition: type = %s',type);
            end
        end

        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary, type)
            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            e_u  = obj.getBoundaryOperator('e',  boundary);
            d1_u = obj.getBoundaryOperator('d1', boundary);
            d2_u = obj.getBoundaryOperator('d2', boundary);
            d3_u = obj.getBoundaryOperator('d3', boundary);
            s_u = obj.getBoundarySign(boundary);

            e_v  = neighbour_scheme.getBoundaryOperator('e',  neighbour_boundary);
            d1_v = neighbour_scheme.getBoundaryOperator('d1', neighbour_boundary);
            d2_v = neighbour_scheme.getBoundaryOperator('d2', neighbour_boundary);
            d3_v = neighbour_scheme.getBoundaryOperator('d3', neighbour_boundary);
            s_v = neighbour_scheme.getBoundarySign(neighbour_boundary);

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

        % Returns the boundary operator op for the boundary specified by the string boundary.
        % op        -- string
        % boundary  -- string
        function o = getBoundaryOperator(obj, op, boundary)
            assertIsMember(op, {'e', 'd1', 'd2', 'd3'})
            assertIsMember(boundary, {'l', 'r'})

            o = obj.([op, '_', boundary]);
        end

        % Returns square boundary quadrature matrix, of dimension
        % corresponding to the number of boundary points
        %
        % boundary -- string
        % Note: for 1d diffOps, the boundary quadrature is the scalar 1.
        function H_b = getBoundaryQuadrature(obj, boundary)
            assertIsMember(boundary, {'l', 'r'})

            H_b = 1;
        end

        % Returns the boundary sign. The right boundary is considered the positive boundary
        % boundary -- string
        function s = getBoundarySign(obj, boundary)
            assertIsMember(boundary, {'l', 'r'})

            switch boundary
                case {'r'}
                    s = 1;
                case {'l'}
                    s = -1;
            end
        end

        function N = size(obj)
            N = obj.grid.N;
        end

    end
end
