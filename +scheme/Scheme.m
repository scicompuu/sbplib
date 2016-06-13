% Start with all matrix returns. When that works see how we should generalize
% to non-matrix stuff/nonlinear
classdef Scheme < handle
    properties (Abstract)
        order % Order accuracy for the approximation

        grid

        D % non-stabalized scheme operator
        H % Discrete norm
    end

    methods (Abstract)
        % Closure functions return the opertors applied to the own doamin to
        % close the boundary Penalty functions return the opertors to force
        % the solution. In the case of an interface it returns the operator
        % applied to the other doamin. In some cases the penalty return value
        % can be ommited and the closure function take care of both parts.
        %       boundary            is a string specifying the boundary e.g.
        %                           'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of
        %                           boundary condition if there are several.
        %       data                is a function returning the data that
        %                           should be applied at the boundary.
        %       neighbour_scheme    is an instance of Scheme that should be
        %                           interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to
        %                           interface to.
        [closure, penalty] = boundary_condition(obj,boundary,type,data)
        [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary)

        % Returns the number of degrees of freedom.
        N = size(obj)
    end

    methods(Static)
        % Calculates the matrcis need for the inteface coupling between
        % boundary bound_u of scheme schm_u and bound_v of scheme schm_v.
        %   [uu, uv, vv, vu] = inteface_coupling(A,'r',B,'l')
        function [uu, uv, vv, vu] = interface_coupling(schm_u,bound_u,schm_v,bound_v)
            [uu,uv] = schm_u.interface(bound_u,schm_v,bound_v);
            [vv,vu] = schm_v.interface(bound_v,schm_u,bound_u);
        end
    end
end
