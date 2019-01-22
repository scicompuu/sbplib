classdef Utux2d < scheme.Scheme
   properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing
        grid % Grid
        order % Order accuracy for the approximation
        v0 % Initial data

        a % Wave speed a = [a1, a2];
          % Can either be a constant vector or a cell array of function handles.

        H % Discrete norm
        H_x, H_y % Norms in the x and y directions
        Hi, Hx, Hy, Hxi, Hyi % Kroneckered norms
        H_w, H_e, H_s, H_n % Boundary quadratures

        % Derivatives
        Dx, Dy

        % Boundary operators
        e_w, e_e, e_s, e_n

        D % Total discrete operator
    end


    methods
         function obj = Utux2d(g ,order, opSet, a)

            default_arg('a',1/sqrt(2)*[1, 1]);
            default_arg('opSet',@sbp.D2Standard);

            assertType(g, 'grid.Cartesian');
            if iscell(a)
                a1 = grid.evalOn(g, a{1});
                a2 = grid.evalOn(g, a{2});
                a = {spdiag(a1), spdiag(a2)};
            else
                a = {a(1), a(2)};
            end

            m = g.size();
            m_x = m(1);
            m_y = m(2);
            m_tot = g.N();

            xlim = {g.x{1}(1), g.x{1}(end)};
            ylim = {g.x{2}(1), g.x{2}(end)};
            obj.grid = g;

            % Operator sets
            ops_x = opSet(m_x, xlim, order);
            ops_y = opSet(m_y, ylim, order);
            Ix = speye(m_x);
            Iy = speye(m_y);

            % Norms
            Hx = ops_x.H;
            Hy = ops_y.H;
            Hxi = ops_x.HI;
            Hyi = ops_y.HI;

            obj.H_w = Hy;
            obj.H_e = Hy;
            obj.H_s = Hx;
            obj.H_n = Hx;
            obj.H_x = Hx;
            obj.H_y = Hy;
            obj.H = kron(Hx,Hy);
            obj.Hi = kron(Hxi,Hyi);
            obj.Hx = kron(Hx,Iy);
            obj.Hy = kron(Ix,Hy);
            obj.Hxi = kron(Hxi,Iy);
            obj.Hyi = kron(Ix,Hyi);

            % Derivatives
            Dx = ops_x.D1;
            Dy = ops_y.D1;
            obj.Dx = kron(Dx,Iy);
            obj.Dy = kron(Ix,Dy);

            % Boundary operators
            obj.e_w = kr(ops_x.e_l, Iy);
            obj.e_e = kr(ops_x.e_r, Iy);
            obj.e_s = kr(Ix, ops_y.e_l);
            obj.e_n = kr(Ix, ops_y.e_r);

            obj.m = m;
            obj.h = [ops_x.h ops_y.h];
            obj.order = order;
            obj.a = a;
            obj.D = -(a{1}*obj.Dx + a{2}*obj.Dy);

        end
        % Closure functions return the opertors applied to the own domain to close the boundary
        % Penalty functions return the opertors to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of boundary condition if there are several.
        %       data                is a function returning the data that should be applied at the boundary.
        %       neighbour_scheme    is an instance of Scheme that should be interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to interface to.
        function [closure, penalty] = boundary_condition(obj,boundary,type)
            default_arg('type','dirichlet');

            sigma = -1; % Scalar penalty parameter
            switch boundary
                case {'w','W','west','West'}
                    tau = sigma*obj.a{1}*obj.e_w*obj.H_y;
                    closure = obj.Hi*tau*obj.e_w';

                case {'s','S','south','South'}
                    tau = sigma*obj.a{2}*obj.e_s*obj.H_x;
                    closure = obj.Hi*tau*obj.e_s';
            end
            penalty = -obj.Hi*tau;

        end

        % type     Struct that specifies the interface coupling.
        %          Fields:
        %          -- couplingType             String, type of interface coupling
        %                                       % Default: 'upwind'. Other: 'centered'
        %          -- interpolation:    type of interpolation, default 'none'
        %          -- interpolationDamping:    damping on upstream and downstream sides, when using interpolation.
        %                                      Default {0,0} gives zero damping.
        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary,type)

            defaultType.couplingType = 'upwind';
            defaultType.interpolation = 'none';
            defaultType.interpolationDamping = {0,0};
            default_struct('type', defaultType);

            switch type.interpolation
            case {'none', ''}
                [closure, penalty] = interfaceStandard(obj,boundary,neighbour_scheme,neighbour_boundary,type);
            case {'op','OP'}
                [closure, penalty] = interfaceNonConforming(obj,boundary,neighbour_scheme,neighbour_boundary,type);
            otherwise
                error('Unknown type of interpolation: %s ', type.interpolation);
            end
        end

        function [closure, penalty] = interfaceStandard(obj,boundary,neighbour_scheme,neighbour_boundary,type)
            couplingType = type.couplingType;

            % Get neighbour boundary operator
            e_neighbour = neighbour_scheme.getBoundaryOperator('e', neighbour_boundary);

            switch couplingType

            % Upwind coupling (energy dissipation)
            case 'upwind'
                 sigma_ds = -1; %"Downstream" penalty
                 sigma_us = 0; %"Upstream" penalty

            % Energy-preserving coupling (no energy dissipation)
            case 'centered'
                 sigma_ds = -1/2; %"Downstream" penalty
                 sigma_us = 1/2; %"Upstream" penalty

            otherwise
                error(['Interface coupling type ' couplingType ' is not available.'])
            end

            switch boundary
                case {'w','W','west','West'}
                    tau = sigma_ds*obj.a{1}*obj.e_w*obj.H_y;
                    closure = obj.Hi*tau*obj.e_w';
                    penalty = -obj.Hi*tau*e_neighbour';
                case {'e','E','east','East'}
                    tau = sigma_us*obj.a{1}*obj.e_e*obj.H_y;
                    closure = obj.Hi*tau*obj.e_e';
                    penalty = -obj.Hi*tau*e_neighbour';
                case {'s','S','south','South'}
                    tau = sigma_ds*obj.a{2}*obj.e_s*obj.H_x;
                    closure = obj.Hi*tau*obj.e_s';
                    penalty = -obj.Hi*tau*e_neighbour';
                case {'n','N','north','North'}
                    tau = sigma_us*obj.a{2}*obj.e_n*obj.H_x;
                    closure = obj.Hi*tau*obj.e_n';
                    penalty = -obj.Hi*tau*e_neighbour';
             end

         end

         function [closure, penalty] = interfaceNonConforming(obj,boundary,neighbour_scheme,neighbour_boundary,type)

            % User can request special interpolation operators by specifying type.interpOpSet
            default_field(type, 'interpOpSet', @sbp.InterpOpsOP);

            interpOpSet = type.interpOpSet;
            couplingType = type.couplingType;
            interpolationDamping = type.interpolationDamping;

            % Get neighbour boundary operator
            e_neighbour = neighbour_scheme.getBoundaryOperator('e', neighbour_boundary);

            switch couplingType

            % Upwind coupling (energy dissipation)
            case 'upwind'
                 sigma_ds = -1; %"Downstream" penalty
                 sigma_us = 0; %"Upstream" penalty

            % Energy-preserving coupling (no energy dissipation)
            case 'centered'
                 sigma_ds = -1/2; %"Downstream" penalty
                 sigma_us = 1/2; %"Upstream" penalty

            otherwise
            error(['Interface coupling type ' couplingType ' is not available.'])
            end

            int_damp_us = interpolationDamping{1};
            int_damp_ds = interpolationDamping{2};

            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            % Find the number of grid points along the interface
            switch boundary
                case {'w','e'}
                    m_u = obj.m(2);
                case {'s','n'}
                    m_u = obj.m(1);
            end
            m_v = size(e_neighbour, 2);

            % Build interpolation operators
            intOps = interpOpSet(m_u, m_v, obj.order, neighbour_scheme.order);
            Iu2v = intOps.Iu2v;
            Iv2u = intOps.Iv2u;

            I_local2neighbour_ds = intOps.Iu2v.bad;
            I_local2neighbour_us = intOps.Iu2v.good;
            I_neighbour2local_ds = intOps.Iv2u.good;
            I_neighbour2local_us = intOps.Iv2u.bad;

            I_back_forth_us = I_neighbour2local_us*I_local2neighbour_us;
            I_back_forth_ds = I_neighbour2local_ds*I_local2neighbour_ds;


            switch boundary
            case {'w','W','west','West'}
                tau = sigma_ds*obj.a{1}*obj.e_w*obj.H_y;
                closure = obj.Hi*tau*obj.e_w';
                penalty = -obj.Hi*tau*I_neighbour2local_ds*e_neighbour';

                beta = int_damp_ds*obj.a{1}...
                        *obj.e_w*obj.H_y;
                closure = closure + obj.Hi*beta*I_back_forth_ds*obj.e_w' - obj.Hi*beta*obj.e_w';
            case {'e','E','east','East'}
                tau = sigma_us*obj.a{1}*obj.e_e*obj.H_y;
                closure = obj.Hi*tau*obj.e_e';
                penalty = -obj.Hi*tau*I_neighbour2local_us*e_neighbour';

                beta = int_damp_us*obj.a{1}...
                        *obj.e_e*obj.H_y;
                closure = closure + obj.Hi*beta*I_back_forth_us*obj.e_e' - obj.Hi*beta*obj.e_e';
            case {'s','S','south','South'}
                tau = sigma_ds*obj.a{2}*obj.e_s*obj.H_x;
                closure = obj.Hi*tau*obj.e_s';
                penalty = -obj.Hi*tau*I_neighbour2local_ds*e_neighbour';

                beta = int_damp_ds*obj.a{2}...
                        *obj.e_s*obj.H_x;
                closure = closure + obj.Hi*beta*I_back_forth_ds*obj.e_s' - obj.Hi*beta*obj.e_s';
            case {'n','N','north','North'}
                tau = sigma_us*obj.a{2}*obj.e_n*obj.H_x;
                closure = obj.Hi*tau*obj.e_n';
                penalty = -obj.Hi*tau*I_neighbour2local_us*e_neighbour';

                beta = int_damp_us*obj.a{2}...
                        *obj.e_n*obj.H_x;
                closure = closure + obj.Hi*beta*I_back_forth_us*obj.e_n' - obj.Hi*beta*obj.e_n';
             end


         end

        % Returns the boundary operator op for the boundary specified by the string boundary.
        % op        -- string
        % boundary  -- string
        function o = getBoundaryOperator(obj, op, boundary)
            assertIsMember(op, {'e'})
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            o = obj.([op, '_', boundary]);
        end

        % Returns square boundary quadrature matrix, of dimension
        % corresponding to the number of boundary points
        %
        % boundary -- string
        function H_b = getBoundaryQuadrature(obj, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            H_b = obj.('H_', boundary);
        end

        function N = size(obj)
            N = obj.m;
        end

    end
end
