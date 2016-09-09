classdef D4Compatible < sbp.OpSet
    properties
        D1 % SBP operator approximating first derivative
        H % Norm matrix
        HI % H^-1
        Q % Skew-symmetric matrix
        e_1 % Left boundary operator
        e_m % Right boundary operator
        D2 % SBP operator for second derivative
        M % Norm matrix, second derivative
        S_1 % Left boundary first derivative
        S_m % Right boundary first derivative
        D3 % SBP operator for third derivative
        Q3 % Skew-symmetric matrix in third derivative
        S2_1 % Left boundary second derivative
        S2_m % Right boundary second derivative
        D4 % SBP operator for fourth derivative
        M4 % Norm matrix, fourth derivative
        S3_1 % Left boundary third derivative
        S3_m % Right boundary third derivative
        m % Number of grid points.
        h % Step size
        x % grid
        borrowing % Struct with borrowing limits for different norm matrices
    end



    methods
        function obj = D4Compatible(m,lim,order)
            
            
            x_l = lim{1};
            x_r = lim{2};
            L = x_r-x_l;
            obj.h = L/(m-1);
            obj.x = linspace(x_l,x_r,m)';

            if order == 2
                [obj.H, obj.HI, obj.D1, obj.D4, obj.e_1, obj.e_m, obj.M4,...
                 obj.Q, obj.S2_1, obj.S2_m, obj.S3_1, obj.S3_m,...
                    obj.S_1, obj.S_m] =...
                    sbp.implementations.d4_compatible_2(m,obj.h);
                obj.borrowing.N.S2 = 0.7500;
                obj.borrowing.N.S3 = 0.3000;
            elseif order == 4
                [obj.H, obj.HI, obj.D1, obj.D4, obj.e_1, obj.e_m, obj.M4,...
                 obj.Q, obj.S2_1, obj.S2_m, obj.S3_1, obj.S3_m,...
                    obj.S_1, obj.S_m] =...
                    sbp.implementations.d4_compatible_4(m,obj.h);
                obj.borrowing.N.S2 = 0.4210;
                obj.borrowing.N.S3 = 0.7080;
            elseif order == 6
                [obj.H, obj.HI, obj.D1, obj.D4, obj.e_1, obj.e_m, obj.M4,...
                 obj.Q, obj.S2_1, obj.S2_m, obj.S3_1, obj.S3_m,...
                    obj.S_1, obj.S_m] =...
                    sbp.implementations.d4_compatible_6(m,obj.h);
                obj.borrowing.N.S2 = 0.06925;
                obj.borrowing.N.S3 = 0.05128;
            else
                error('Invalid operator order.');
            end

            obj.m = m;

            obj.D2 = [];
            obj.D3 = [];


        end
    end



end