classdef Utux < scheme.Scheme
   properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing
        grid % Grid
        order % Order accuracy for the approximation

        H % Discrete norm
        D

        D1
        Hi
        e_l
        e_r
        v0
    end


    methods
        function obj = Utux(g, order, opSet)
            default_arg('opSet',@sbp.D2Standard);

            m = g.size();
            xl = g.getBoundary('l');
            xr = g.getBoundary('r');
            xlim = {xl, xr};

            ops = opSet(m, xlim, order);
            obj.D1 = ops.D1;

            obj.grid = g;

            obj.H =  ops.H;
            obj.Hi = ops.HI;

            obj.e_l = ops.e_l;
            obj.e_r = ops.e_r;
            obj.D = -obj.D1;

            obj.m = m;
            obj.h = ops.h;
            obj.order = order;

        end
        % Closure functions return the opertors applied to the own doamin to close the boundary
        % Penalty functions return the opertors to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of boundary condition if there are several.
        %       data                is a function returning the data that should be applied at the boundary.
        %       neighbour_scheme    is an instance of Scheme that should be interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to interface to.
        function [closure, penalty] = boundary_condition(obj,boundary,type)
            default_arg('type','dirichlet');
            tau =-1*obj.e_l;
            closure = obj.Hi*tau*obj.e_l';
            penalty = -obj.Hi*tau;

         end

         function [closure, penalty] = interface(obj, boundary, neighbour_scheme, neighbour_boundary, type)
             switch boundary
                 % Upwind coupling
                 case {'l','left'}
                     tau = -1*obj.e_l;
                     closure = obj.Hi*tau*obj.e_l';
                     penalty = -obj.Hi*tau*neighbour_scheme.e_r';
                 case {'r','right'}
                     tau = 0*obj.e_r;
                     closure = obj.Hi*tau*obj.e_r';
                     penalty = -obj.Hi*tau*neighbour_scheme.e_l';
             end

         end

        % Returns the boundary operator op for the boundary specified by the string boundary.
        % op        -- string or a cell array of strings
        % boundary  -- string
        function varargout = getBoundaryOperator(obj, op, boundary)

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
                    otherwise
                        error('No such boundary: boundary = %s',boundary);
                    end
                    varargout{i} = e;
                end
            end
        end

        function N = size(obj)
            N = obj.m;
        end

    end

    methods(Static)
        % Calculates the matrices needed for the inteface coupling between boundary bound_u of scheme schm_u
        % and bound_v of scheme schm_v.
        %   [uu, uv, vv, vu] = inteface_coupling(A,'r',B,'l')
        function [uu, uv, vv, vu] = interface_coupling(schm_u,bound_u,schm_v,bound_v)
            [uu,uv] = schm_u.interface(bound_u,schm_v,bound_v);
            [vv,vu] = schm_v.interface(bound_v,schm_u,bound_u);
        end
    end
end