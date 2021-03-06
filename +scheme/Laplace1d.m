classdef Laplace1d < scheme.Scheme
    properties
        grid
        order % Order accuracy for the approximation

        D % non-stabalized scheme operator
        H % Discrete norm
        M % Derivative norm
        a

        D2
        Hi
        e_l
        e_r
        d_l
        d_r
        gamm
    end

    methods
        function obj = Laplace1d(grid, order, a)
            default_arg('a', 1);

            assertType(grid, 'grid.Cartesian');

            ops = sbp.D2Standard(grid.size(), grid.lim{1}, order);

            obj.D2 = sparse(ops.D2);
            obj.H =  sparse(ops.H);
            obj.Hi = sparse(ops.HI);
            obj.M =  sparse(ops.M);
            obj.e_l = sparse(ops.e_l);
            obj.e_r = sparse(ops.e_r);
            obj.d_l = -sparse(ops.d1_l);
            obj.d_r = sparse(ops.d1_r);


            obj.grid = grid;
            obj.order = order;

            obj.a = a;
            obj.D = a*obj.D2;

            obj.gamm = grid.h*ops.borrowing.M.S;
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

            e = obj.getBoundaryOperator('e', boundary);
            d = obj.getBoundaryOperator('d', boundary);
            s = obj.getBoundarySign(boundary);

            switch type
                % Dirichlet boundary condition
                case {'D','d','dirichlet'}
                    tuning = 1.1;
                    tau1 = -tuning/obj.gamm;
                    tau2 =  1;

                    tau = tau1*e + tau2*d;

                    closure = obj.a*obj.Hi*tau*e';
                    penalty = -obj.a*obj.Hi*tau;

                % Neumann boundary condition
                case {'N','n','neumann'}
                    tau = -e;

                    closure = obj.a*obj.Hi*tau*d';
                    penalty = -obj.a*obj.Hi*tau;

                % Unknown, boundary condition
                otherwise
                    error('No such boundary condition: type = %s',type);
            end
        end

        function [closure, penalty] = interface(obj, boundary, neighbour_scheme, neighbour_boundary, type)
            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            e_u = obj.getBoundaryOperator('e', boundary);
            d_u = obj.getBoundaryOperator('d', boundary);
            s_u = obj.getBoundarySign(boundary);

            e_v = neighbour_scheme.getBoundaryOperator('e', neighbour_boundary);
            d_v = neighbour_scheme.getBoundaryOperator('d', neighbour_boundary);
            s_v = neighbour_scheme.getBoundarySign(neighbour_boundary);

            a_u = obj.a;
            a_v = neighbour_scheme.a;

            gamm_u = obj.gamm;
            gamm_v = neighbour_scheme.gamm;

            tuning = 1.1;

            tau1 = -1/4*(a_u/gamm_u + a_v/gamm_v) * tuning;
            tau2 = 1/2*a_u;
            sig1 = -1/2;
            sig2 = 0;

            tau = tau1*e_u + tau2*d_u;
            sig = sig1*e_u + sig2*d_u;

            closure = obj.Hi*( tau*e_u' + sig*a_u*d_u');
            penalty = obj.Hi*(-tau*e_v' + sig*a_v*d_v');
        end

        % Returns the boundary operator op for the boundary specified by the string boundary.
        % op        -- string
        % boundary  -- string
        function o = getBoundaryOperator(obj, op, boundary)
            assertIsMember(op, {'e', 'd'})
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
            N = obj.grid.size();
        end

    end
end
