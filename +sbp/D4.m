classdef D4 < sbp.OpSet
    properties
        norms % Struct containing norm matrices such as H,Q, M
        boundary  % Struct contanging vectors for boundry point approximations
        derivatives % Struct containging differentiation operators
        borrowing % Struct with borrowing limits for different norm matrices
        m % Number of grid points.
        h % Step size
    end



    methods
        function obj = D4(m,h,order)

            if order == 4
                [H, HI, D1, D2, D3, D4, e_1, e_m, M, M4,Q, Q3, S2_1,...
                    S2_m, S3_1, S3_m, S_1, S_m] = sbp.implementations.d4_4(m,h);
                obj.borrowing.N.S2 = 0.5485;
                obj.borrowing.N.S3 = 1.0882;
            elseif order == 6
                [H, HI, D1, D2, D3, D4, e_1, e_m, M, M4,Q, Q3, S2_1,...
                    S2_m, S3_1, S3_m, S_1, S_m] = sbp.implementations.d4_6(m,h);
                obj.borrowing.N.S2 = 0.3227;
                obj.borrowing.N.S3 = 0.1568;
            else
                error('Invalid operator order %d.',order);
            end

            obj.h = h;
            obj.m = m;

            obj.norms.H = H;
            obj.norms.HI = HI;
            obj.norms.Q = Q;
            obj.norms.M = M;
            obj.norms.Q3 = Q3;
            obj.norms.N = M4;

            obj.boundary.e_1 = e_1;
            obj.boundary.S_1 = S_1;
            obj.boundary.S2_1 = S2_1;
            obj.boundary.S3_1 = S3_1;

            obj.boundary.e_m = e_m;
            obj.boundary.S_m = S_m;
            obj.boundary.S2_m = S2_m;
            obj.boundary.S3_m = S3_m;

            obj.derivatives.D1 = D1;
            obj.derivatives.D2 = D2;
            obj.derivatives.D3 = D3;
            obj.derivatives.D4 = D4;

        end
    end

    methods (Static)
        function lambda = smallestGrid(obj)
            error('Not implmented')
        end
    end



end