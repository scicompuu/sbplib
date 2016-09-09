classdef D2BlockNorm < sbp.OpSet
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
        m % Number of grid points.
        h % Step size
        x % grid
        borrowing % Struct with borrowing limits for different norm matrices
    end



    methods
        function obj = D2BlockNorm(m,lim,order)

            x_l = lim{1};
            x_r = lim{2};
            L = x_r-x_l;
            obj.h = L/(m-1);
            obj.x = linspace(x_l,x_r,m)';

            if order == 4
                [obj.H, obj.HI, obj.D1, obj.D2, obj.e_1,...
                 obj.e_m, obj.M, obj.Q, obj.S_1, obj.S_m] = ...
                    sbp.implementations.d2_blocknorm_4(m,obj.h);
            elseif order == 6
                [obj.H, obj.HI, obj.D1, obj.D2, obj.e_1,...
                 obj.e_m, obj.M, obj.Q, obj.S_1, obj.S_m] = ...
                    sbp.implementations.d2_blocknorm_6(m,obj.h);
            elseif order == 8
                [obj.H, obj.HI, obj.D1, obj.D2, obj.e_1,...
                 obj.e_m, obj.M, obj.Q, obj.S_1, obj.S_m] = ...
                    sbp.implementations.d2_blocknorm_8(m,obj.h);
            elseif order == 10
                [obj.H, obj.HI, obj.D1, obj.D2, obj.e_1,...
                 obj.e_m, obj.M, obj.Q, obj.S_1, obj.S_m] = ...
                    sbp.implementations.d2_blocknorm_10(m,obj.h);
            else
                error('Invalid operator order %d.',order);
            end

            obj.m = m;

        end
    end



end