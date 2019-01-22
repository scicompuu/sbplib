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

            [e, d] = obj.getBoundaryOperator({'e', 'd'}, boundary);
            s = obj.getBoundarySign(boundary);

            switch type
                % Dirichlet boundary condition
                case {'D','dirichlet'}
                    tuning = 1.1;
                    tau1 = -tuning/obj.gamm;
                    tau2 =  1;

                    tau = tau1*e + tau2*d;

                    closure = obj.a*obj.Hi*tau*e';
                    penalty = obj.a*obj.Hi*tau;

                % Neumann boundary condition
                case {'N','neumann'}
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
            [e_u, d_u] = obj.getBoundaryOperator({'e', 'd'}, boundary);
            s_u = obj.getBoundarySign(boundary);

            [e_v, d_v] = neighbour_scheme.getBoundaryOperator({'e', 'd'}, neighbour_boundary);
            s_v = neighbour_scheme.getBoundarySign(neighbour_boundary);

            a_u = obj.a;
            a_v = neighbour_scheme.a;

            gamm_u = obj.gamm;
            gamm_v = neighbour_scheme.gamm;

            tuning = 1.1;

            tau1 = -(a_u/gamm_u + a_v/gamm_v) * tuning;
            tau2 = 1/2*a_u;
            sig1 = -1/2;
            sig2 = 0;

            tau = tau1*e_u + tau2*d_u;
            sig = sig1*e_u + sig2*d_u;

            closure = obj.Hi*( tau*e_u' + sig*a_u*d_u');
            penalty = obj.Hi*(-tau*e_v' + sig*a_v*d_v');
        end

        % Returns the boundary operator op for the boundary specified by the string boundary.
        % op        -- string or a cell array of strings
        % boundary  -- string
        function varargout = getBoundaryOperator(obj, op, boundary)
            assertIsMember(boundary, {'l', 'r'})

            if ~iscell(op)
                op = {op};
            end

            for i = 1:numel(op)
                switch op{i}
                case 'e'
                    switch boundary
                    case 'l'
                        e = obj.e_l;
                    case 'r'
                        e = obj.e_r;
                    end
                    varargout{i} = e;

                case 'd'
                    switch boundary
                    case 'l'
                        d = obj.d_l;
                    case 'r'
                        d = obj.d_r;
                    end
                    varargout{i} = d;
                end
            end
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