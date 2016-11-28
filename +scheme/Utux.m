classdef Utux < scheme.Scheme
   properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing
        x % Grid
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
         function obj = Utux(m,xlim,order)
             default_arg('a',1);
           
           %Old operators  
           % [x, h] = util.get_grid(xlim{:},m);
           %ops = sbp.Ordinary(m,h,order);
           
             % ops = sbp.D1Nonequidistant(m,xlim,order);
              ops = sbp.D2Standard(m,xlim,order);
             obj.D1 = ops.D1;
%              ops = sbp.D1Upwind(m,xlim,order);
%             obj.D1 = ops.Dm;
            obj.x=ops.x;

            
            obj.H =  ops.H;
            obj.Hi = ops.HI;
        
            obj.e_l = ops.e_l;
            obj.e_r = ops.e_r;
            obj.D=obj.D1;

            obj.m = m;
            obj.h = ops.h;
            obj.order = order;
            obj.x = ops.x;

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
            tau =-1*obj.e_l;  
            closure = obj.Hi*tau*obj.e_l';       
            penalty = 0*obj.e_l;
                
         end
          
         function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary)
          error('An interface function does not exist yet');
         end
      
        function N = size(obj)
            N = obj.m;
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