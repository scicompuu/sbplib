classdef D2 < sbp.OpSet
    properties
        norms % Struct containing norm matrices such as H,Q, M
        boundary  % Struct contanging vectors for boundry point approximations
        derivatives % Struct containging differentiation operators
        borrowing % Struct with borrowing limits for different norm matrices
        m % Number of grid points.
        h % Step size
    end

    methods
        function obj = D2(m,h,order)

            if order == 2
                [H, HI, D1, D2, e_1, e_m, M,Q S_1, S_m] = sbp.ordinary2(m,h);
                obj.borrowing.M.S = 0.4000;
            elseif order == 4
                [H, HI, D1, D2, e_1, e_m, M,Q S_1, S_m] = sbp.ordinary4(m,h);
                obj.borrowing.M.S = 0.2508;
            elseif order == 6
                [H, HI, D1, D2, e_1, e_m, M,Q S_1, S_m] = sbp.ordinary6(m,h);
                obj.borrowing.M.S = 0.1878;
            elseif order == 8
                [H, HI, D1, D2, e_1, e_m, M,Q S_1, S_m] = sbp.ordinary8(m,h);
                obj.borrowing.M.S = 0.0015;
            elseif order == 10
                [H, HI, D1, D2, e_1, e_m, M,Q S_1, S_m] = sbp.ordinary10(m,h);
                obj.borrowing.M.S = 0.0351;
            else
                error('Invalid operator order %d.',order);
            end

            obj.h = h;
            obj.m = m;

            obj.norms.H = H;
            obj.norms.HI = HI;
            obj.norms.Q = Q;
            obj.norms.M = M;

            obj.boundary.e_1 = e_1;
            obj.boundary.S_1 = S_1;

            obj.boundary.e_m = e_m;
            obj.boundary.S_m = S_m;

            obj.derivatives.D1 = D1;
            obj.derivatives.D2 = D2;

        end
    end

    methods (Static)
        function lambda = smallestGrid(obj)
            error('Not implmented')
        end
    end
end





