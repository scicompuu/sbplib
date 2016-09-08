classdef D4CompatibleVariable < sbp.OpSet
    properties
        norms % Struct containing norm matrices such as H,Q, M
        boundary  % Struct contanging vectors for boundry point approximations
        derivatives % Struct containging differentiation operators
        borrowing % Struct with borrowing limits for different norm matrices
        m % Number of grid points.
        h % Step size
    end



    methods
        function obj = D4CompatibleVariable(m,h,order)

            if order == 2
                [H, HI, ~, D2, ~, D4, e_1, e_m, M4, ~, S2_1, S2_m, S3_1,...
                    S3_m, S_1, S_m] =...
                    sbp.implementations.d4_compatible_halfvariable_2(m,h);
                obj.borrowing.N.S2 = 1.2500;
                obj.borrowing.N.S3 = 0.4000;
            elseif order == 4
                [H, HI, D2, D4, e_1, e_m, M4, S2_1, S2_m, S3_1, S3_m, S_1,...
                    S_m] =...
                    sbp.implementations.d4_compatible_halfvariable_4(m,h);
                obj.borrowing.N.S2 = 0.5055;
                obj.borrowing.N.S3 = 0.9290;
            elseif order == 6
                [H, HI, D2, D4, e_1, e_m, M4, S2_1, S2_m, S3_1, S3_m, S_1,...
                    S_m] =...
                    sbp.implementations.d4_compatible_halfvariable_6(m,h);
                obj.borrowing.N.S2 = 0.3259;
                obj.borrowing.N.S3 = 0.1580;
            else
                error('Invalid operator order.');
            end

            obj.h = h;
            obj.m = m;

            obj.norms.H = H;
            obj.norms.HI = HI;
            obj.norms.N = M4;

            obj.boundary.e_1 = e_1;
            obj.boundary.S_1 = S_1;
            obj.boundary.S2_1 = S2_1;
            obj.boundary.S3_1 = S3_1;

            obj.boundary.e_m = e_m;
            obj.boundary.S_m = S_m;
            obj.boundary.S2_m = S2_m;
            obj.boundary.S3_m = S3_m;

            obj.derivatives.D2 = D2;
            obj.derivatives.D4 = D4;

        end
    end



end