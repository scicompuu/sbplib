classdef D1Upwind < sbp.OpSet
    properties
        norms % Struct containing norm matrices such as H,Q, M
        boundary  % Struct contanging vectors for boundry point approximations
        derivatives % Struct containging differentiation operators
        borrowing % Struct with borrowing limits for different norm matrices
        m % Number of grid points.
        h % Step size
    end

    methods
        function obj = D1Upwind(m,h,order)

            switch order
                case 2
                    [H, HI, Dp, Dm, e_1, e_m] = sbp.upwind2(m,h);
                case 3
                    [H, HI, Dp, Dm, e_1, e_m] = sbp.upwind3(m,h);
                case 4
                    [H, HI, Dp, Dm, e_1, e_m] = sbp.upwind4(m,h);
                case 5
                    [H, HI, Dp, Dm, e_1, e_m] = sbp.upwind5(m,h);
                case 6
                    [H, HI, Dp, Dm, e_1, e_m] = sbp.upwind6(m,h);
                case 7
                    [H, HI, Dp, Dm, e_1, e_m] = sbp.upwind7(m,h);
                case 8
                    [H, HI, Dp, Dm, e_1, e_m] = sbp.upwind8(m,h);
                case 9
                    [H, HI, Dp, Dm, e_1, e_m] = sbp.upwind9(m,h);
                otherwise
                    error('Invalid operator order %d.',order);
            end

            obj.h = h;
            obj.m = m;

            obj.norms.H = H;
            obj.norms.HI = HI;

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





