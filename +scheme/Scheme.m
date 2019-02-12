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
        %       neighbour_scheme    is an instance of Scheme that should be
        %                           interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to
        %                           interface to.
        %       penalty  may be a cell array if there are several penalties with different weights
        [closure, penalty] = boundary_condition(obj,boundary,type) % TODO: Change name to boundaryCondition

        % type -- sets the type of interface, could be a string or a struct or something else
        %         depending on the particular scheme implementation
        [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary,type)

        op = getBoundaryOperator(obj, opName, boundary)
        H_b= getBoundaryQuadrature(obj, boundary)

        % Returns the number of degrees of freedom.
        N = size(obj)
    end
end
