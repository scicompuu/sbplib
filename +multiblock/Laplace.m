classdef Laplace < scheme.Scheme
    properties
        grid
        order
        mbDiffOp

        D
        H
        J
    end
    methods
        function obj = Laplace(g, order, a, b, opGen)
            default_arg('order', 4);
            default_arg('a', 1);
            default_arg('b', 1);
            default_arg('opGen', @sbp.D4Variable);

            obj.grid = g;
            obj.order = order;
            obj.mbDiffOp = multiblock.DiffOp(@scheme.LaplaceCurvilinear, obj.grid, order, {a,b,opGen});

            obj.D = obj.mbDiffOp.D;
            obj.J = obj.jacobian();
            obj.H = obj.mbDiffOp.H * obj.jacobian();
        end

        function s = size(obj)
            s = size(obj.mbDiffOp);
        end

        function J = jacobian(obj)
            N = obj.grid.nBlocks;
            J = cell(N,N);

            for i = 1:N
                J{i,i} = obj.mbDiffOp.diffOps{i}.J;
            end
            J = blockmatrix.toMatrix(J);
        end

        function op = getBoundaryOperator(obj, opName, boundary)
            op = getBoundaryOperator(obj.mbDiffOp, opName, boundary);
        end

        function op = getBoundaryQuadrature(obj, boundary)
            op = getBoundaryQuadrature(obj.mbDiffOp, boundary);
        end

        function [closure, penalty] = boundary_condition(obj,boundary,type) % TODO: Change name to boundaryCondition
            [closure, penalty] = boundary_condition(obj.mbDiffOp, boundary, type);
        end
        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary)
            error('Not implemented')
        end
    end
end
