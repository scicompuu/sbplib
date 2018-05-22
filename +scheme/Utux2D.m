classdef Utux2D < scheme.Scheme
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

        % Derivatives
        Dx, Dy
        
        % Boundary operators
        e_w, e_e, e_s, e_n
        
        D % Total discrete operator

        % String, type of interface coupling
        % Default: 'upwind'
        % Other:   'centered'
        coupling_type 

        % String, type of interpolation operators
        % Default: 'AWW' (Almquist Wang Werpers)
        % Other:   'MC' (Mattsson Carpenter)
        interpolation_type

        
        % Cell array, damping on upwstream and downstream sides.
        interpolation_damping

    end


    methods 
         function obj = Utux2D(g ,order, opSet, a, coupling_type, interpolation_type, interpolation_damping)
            
            default_arg('interpolation_damping',{0,0});
            default_arg('interpolation_type','AWW'); 
            default_arg('coupling_type','upwind'); 
            default_arg('a',1/sqrt(2)*[1, 1]); 
            default_arg('opSet',@sbp.D2Standard);

            assert(isa(g, 'grid.Cartesian'))
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
            obj.interpolation_type = interpolation_type;
            obj.interpolation_damping = interpolation_damping;
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
          
         function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary)
             
             % Get neighbour boundary operator
             switch neighbour_boundary
                 case {'e','E','east','East'}
                     e_neighbour = neighbour_scheme.e_e;
                     m_neighbour = neighbour_scheme.m(2);
                 case {'w','W','west','West'}
                     e_neighbour = neighbour_scheme.e_w;
                     m_neighbour = neighbour_scheme.m(2);
                 case {'n','N','north','North'}
                     e_neighbour = neighbour_scheme.e_n;
                     m_neighbour = neighbour_scheme.m(1);
                 case {'s','S','south','South'}
                     e_neighbour = neighbour_scheme.e_s;
                     m_neighbour = neighbour_scheme.m(1);
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

             % Check grid ratio for interpolation
             switch boundary
                 case {'w','W','west','West','e','E','east','East'}
                     m = obj.m(2);       
                 case {'s','S','south','South','n','N','north','North'}
                     m = obj.m(1);
             end
             grid_ratio = m/m_neighbour;
             if grid_ratio ~= 1

                [ms, index] = sort([m, m_neighbour]);
                orders = [obj.order, neighbour_scheme.order];
                orders = orders(index);

                switch obj.interpolation_type
                case 'MC'
                    interpOpSet = sbp.InterpMC(ms(1),ms(2),orders(1),orders(2));
                    if grid_ratio < 1
                        I_neighbour2local_us = interpOpSet.IF2C;
                        I_neighbour2local_ds = interpOpSet.IF2C;
                        I_local2neighbour_us = interpOpSet.IC2F;
                        I_local2neighbour_ds = interpOpSet.IC2F;
                    elseif grid_ratio > 1
                        I_neighbour2local_us = interpOpSet.IC2F;
                        I_neighbour2local_ds = interpOpSet.IC2F;
                        I_local2neighbour_us = interpOpSet.IF2C;
                        I_local2neighbour_ds = interpOpSet.IF2C;
                    end
                case 'AWW'
                    %String 'C2F' indicates that ICF2 is more accurate.
                    interpOpSetF2C = sbp.InterpAWW(ms(1),ms(2),orders(1),orders(2),'F2C');
                    interpOpSetC2F = sbp.InterpAWW(ms(1),ms(2),orders(1),orders(2),'C2F'); 
                    if grid_ratio < 1 
                        % Local is coarser than neighbour
                        I_neighbour2local_us = interpOpSetC2F.IF2C;
                        I_neighbour2local_ds = interpOpSetF2C.IF2C;
                        I_local2neighbour_us = interpOpSetC2F.IC2F;
                        I_local2neighbour_ds = interpOpSetF2C.IC2F;
                    elseif grid_ratio > 1
                        % Local is finer than neighbour 
                        I_neighbour2local_us = interpOpSetF2C.IC2F;
                        I_neighbour2local_ds = interpOpSetC2F.IC2F;
                        I_local2neighbour_us = interpOpSetF2C.IF2C;
                        I_local2neighbour_ds = interpOpSetC2F.IF2C;
                    end
                otherwise
                    error(['Interpolation type ' obj.interpolation_type ...
                         ' is not available.' ]);
                end

             else 
                % No interpolation required
                I_neighbour2local_us = speye(m,m);
                I_neighbour2local_ds = speye(m,m);
            end    
             
             int_damp_us = obj.interpolation_damping{1};
             int_damp_ds = obj.interpolation_damping{2};

             I = speye(m,m);
             I_back_forth_us = I_neighbour2local_us*I_local2neighbour_us;
             I_back_forth_ds = I_neighbour2local_ds*I_local2neighbour_ds;


             switch boundary
                 case {'w','W','west','West'}
                     tau = sigma_ds*obj.a{1}*obj.e_w*obj.H_y;
                     closure = obj.Hi*tau*obj.e_w';
                     penalty = -obj.Hi*tau*I_neighbour2local_ds*e_neighbour';

                     beta = int_damp_ds*obj.a{1}...
                            *obj.e_w*obj.H_y;
                     closure = closure + obj.Hi*beta*(I_back_forth_ds - I)*obj.e_w';     
                 case {'e','E','east','East'}
                     tau = sigma_us*obj.a{1}*obj.e_e*obj.H_y;
                     closure = obj.Hi*tau*obj.e_e';
                     penalty = -obj.Hi*tau*I_neighbour2local_us*e_neighbour';

                     beta = int_damp_us*obj.a{1}...
                            *obj.e_e*obj.H_y;
                     closure = closure + obj.Hi*beta*(I_back_forth_us - I)*obj.e_e'; 
                 case {'s','S','south','South'}
                     tau = sigma_ds*obj.a{2}*obj.e_s*obj.H_x;
                     closure = obj.Hi*tau*obj.e_s'; 
                     penalty = -obj.Hi*tau*I_neighbour2local_ds*e_neighbour';

                     beta = int_damp_ds*obj.a{2}...
                            *obj.e_s*obj.H_x;
                     closure = closure + obj.Hi*beta*(I_back_forth_ds - I)*obj.e_s';
                 case {'n','N','north','North'}
                     tau = sigma_us*obj.a{2}*obj.e_n*obj.H_x;
                     closure = obj.Hi*tau*obj.e_n';
                     penalty = -obj.Hi*tau*I_neighbour2local_us*e_neighbour';

                     beta = int_damp_us*obj.a{2}...
                            *obj.e_n*obj.H_x;
                     closure = closure + obj.Hi*beta*(I_back_forth_us - I)*obj.e_n'; 
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