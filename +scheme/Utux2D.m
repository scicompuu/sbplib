classdef Utux2D < scheme.Scheme
   properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing
        grid % Grid
        order % Order accuracy for the approximation
        v0 % Initial data
        
        a % Wave speed a = [a1, a2];

        H % Discrete norm
        H_x, H_y % Norms in the x and y directions
        Hi, Hx, Hy, Hxi, Hyi % Kroneckered norms

        % Derivatives
        Dx, Dy
        
        % Boundary operators
        e_w, e_e, e_s, e_n
        
        D % Total discrete operator

        % String, type of interface coupling
        % Default: 'upwind'
        % Other:   'centered'
        coupling_type 

        
    end


    methods 
         function obj = Utux2D(g ,order, opSet, a, coupling_type)
            
            default_arg('coupling_type','upwind'); 
            default_arg('a',1/sqrt(2)*[1, 1]); 
            default_arg('opSet',@sbp.D2Standard);
            assert(isa(g, 'grid.Cartesian'))
             
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
            obj.coupling_type = coupling_type;
            obj.D = -(a(1)*obj.Dx + a(2)*obj.Dy);

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
                    tau = sigma*obj.a(1)*obj.e_w*obj.H_y;
                    closure = obj.Hi*tau*obj.e_w';
                    
                case {'s','S','south','South'}
                    tau = sigma*obj.a(2)*obj.e_s*obj.H_x;
                    closure = obj.Hi*tau*obj.e_s';
            end  
            penalty = -obj.Hi*tau;
                
         end
          
         function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary)
             
             % Get neighbour boundary operator
             switch neighbour_boundary
                 case {'e','E','east','East'}
                     e_neighbour = neighbour_scheme.e_e;
                 case {'w','W','west','West'}
                     e_neighbour = neighbour_scheme.e_w;
                 case {'n','N','north','North'}
                     e_neighbour = neighbour_scheme.e_n;
                 case {'s','S','south','South'}
                     e_neighbour = neighbour_scheme.e_s;
             end
             
             switch obj.coupling_type
             
             % Upwind coupling (energy dissipation)
             case 'upwind'
                 sigma_ds = -1; %"Downstream" penalty
                 sigma_us = 0; %"Upstream" penalty

             % Energy-preserving coupling (no energy dissipation)
             case 'centered'
                 sigma_ds = -1/2; %"Downstream" penalty
                 sigma_us = 1/2; %"Upstream" penalty

             otherwise
                error(['Interface coupling type ' coupling_type ' is not available.'])
             end
             
             switch boundary
                 case {'w','W','west','West'}
                     tau = sigma_ds*obj.a(1)*obj.e_w*obj.H_y;
                     closure = obj.Hi*tau*obj.e_w';       
                 case {'e','E','east','East'}
                     tau = sigma_us*obj.a(1)*obj.e_e*obj.H_y;
                     closure = obj.Hi*tau*obj.e_e';
                 case {'s','S','south','South'}
                     tau = sigma_ds*obj.a(2)*obj.e_s*obj.H_x;
                     closure = obj.Hi*tau*obj.e_s'; 
                 case {'n','N','north','North'}
                     tau = sigma_us*obj.a(2)*obj.e_n*obj.H_x;
                     closure = obj.Hi*tau*obj.e_n';
             end
             penalty = -obj.Hi*tau*e_neighbour';
                 
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