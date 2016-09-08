classdef D2Variable < sbp.OpSet
    properties
        norms % Struct containing norm matrices such as H,Q, M
        boundary  % Struct contanging vectors for boundry point approximations
        derivatives % Struct containging differentiation operators
        borrowing % Struct with borrowing limits for different norm matrices
        m % Number of grid points.
        h % Step size
    end

    methods
        function obj = D2Variable(m,h,order)

            switch order
                case 4
                    [H, HI, D1, D2, e_1, e_m, S_1, S_m] = ...
                        sbp.implementations.d2_variable_4(m,h);
                    obj.borrowing.M.S = 0.2505765857;
                otherwise
                    error('Invalid operator order %d.',order);
            end

            obj.h = h;
            obj.m = m;

            obj.norms.H = H;
            obj.norms.HI = HI;
            % obj.norms.Q = Q;
            % obj.norms.M = M;

            obj.boundary.e_1 = e_1;
            obj.boundary.S_1 = S_1;

            obj.boundary.e_m = e_m;
            obj.boundary.S_m = S_m;

            obj.derivatives.D1 = D1;
            obj.derivatives.D2 = D2;

        end
    end

end





