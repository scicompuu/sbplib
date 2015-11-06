classdef Ordinary < sbp.OpSet
    properties
        norms % Struct containing norm matrices such as H,Q, M
        boundary  % Struct contanging vectors for boundry point approximations
        derivatives % Struct containging differentiation operators
        borrowing % Struct with borrowing limits for different norm matrices
        m % Number of grid points.
        h % Step size
    end

    methods
        function obj = Ordinary(m,h,order)

            if order == 3
                [H, HI, Dp, Dm, e_1, e_m] = sbp.upwind3(m,h);
            elseif order == 4
                [H, HI, Dp, Dm, e_1, e_m] = sbp.upwind4(m,h);
            elseif order == 6
                [H, HI, Dp, Dm, e_1, e_m] = sbp.upwind6(m,h);
            elseif order == 8
                [H, HI, Dp, Dm, e_1, e_m] = sbp.upwind8(m,h);
            else
                error('Invalid operator order %d.',order);
            end

            obj.h = h;
            obj.m = m;

            obj.norms.H = H;
            obj.norms.HI = HI;
            obj.norms.Q = Q;

            obj.boundary.e_1 = e_1;
            obj.boundary.e_m = e_m;

            obj.derivatives.Dp = Dp;
            obj.derivatives.Dm = Dm;
        end
    end

    methods (Static)
        function lambda = smallestGrid(obj)
            error('Not implmented')
        end
    end
end





