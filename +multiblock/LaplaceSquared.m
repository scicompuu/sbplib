classdef LaplaceSquared < scheme.Scheme
    properties
        grid
        order
        laplaceDiffOp

        D
        H
        Hi

        a,b
    end

    methods
        % Discretisation of a*nabla*b*nabla
        function obj = LaplaceSquared(g, order, a, b, opGen)
            default_arg('order', 4);
            default_arg('a', 1);
            default_arg('b', 1);
            default_arg('opGen', @sbp.D4Variable);

            if isscalar(a)
                a = grid.evalOn(g, a);
            end

            if isscalar(b)
                b = grid.evalOn(g, b);
            end

            obj.grid = g;
            obj.order = order;
            obj.a = a;
            obj.b = b;

            obj.laplaceDiffOp = multiblock.Laplace(g, order, 1, 1, opGen);

            obj.H = obj.laplaceDiffOp.H;
            obj.Hi = spdiag(1./diag(obj.H));

            A = spdiag(a);
            B = spdiag(b);

            D_laplace = obj.laplaceDiffOp.D;
            obj.D = A*D_laplace*B*D_laplace;
        end

        function s = size(obj)
            s = size(obj.laplaceDiffOp);
        end

        function op = getBoundaryOperator(obj, opName, boundary)
            switch opName
                case 'e'
                    op = getBoundaryOperator(obj.laplaceDiffOp, 'e', boundary);
                case 'd1'
                    op = getBoundaryOperator(obj.laplaceDiffOp, 'd', boundary);
                case 'd2'
                    e = getBoundaryOperator(obj.laplaceDiffOp, 'e', boundary);
                    op = (e'*obj.laplaceDiffOp.D)';
                case 'd3'
                    d1 = getBoundaryOperator(obj.laplaceDiffOp, 'd', boundary);
                    op = (d1'*spdiag(obj.b)*obj.laplaceDiffOp.D)';
            end
        end

        function op = getBoundaryQuadrature(obj, boundary)
            op = getBoundaryQuadrature(obj.laplaceDiffOp, boundary);
        end

        function [closure, penalty] = boundary_condition(obj,boundary,type) % TODO: Change name to boundaryCondition
            switch type
                case 'e'
                    error('Bc of type ''e'' not implemented')
                case 'd1'
                    error('Bc of type ''d1'' not implemented')
                case 'd2'
                    e = obj.getBoundaryOperator('e', boundary);
                    d1 = obj.getBoundaryOperator('d1', boundary);
                    d2 = obj.getBoundaryOperator('d2', boundary);
                    H_b = obj.getBoundaryQuadrature(boundary);

                    A = spdiag(obj.a);
                    B_b = spdiag(e'*obj.b);

                    tau = obj.Hi*A*d1*B_b*H_b;
                    closure =  tau*d2';
                    penalty = -tau;
                case 'd3'
                    e = obj.getBoundaryOperator('e', boundary);
                    d3 = obj.getBoundaryOperator('d1', boundary);
                    H_b = obj.getBoundaryQuadrature(boundary);

                    A = spdiag(obj.a);

                    tau = -obj.Hi*A*e*H_b;
                    closure =  tau*d3';
                    penalty = -tau;
            end
        end

        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary)
            error('Not implemented')
        end
    end
end