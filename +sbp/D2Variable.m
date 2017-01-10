classdef D2Variable < sbp.OpSet
    properties
        D1 % SBP operator approximating first derivative
        H % Norm matrix
        HI % H^-1
        Q % Skew-symmetric matrix
        e_l % Left boundary operator
        e_r % Right boundary operator
        D2 % SBP operator for second derivative
        M % Norm matrix, second derivative
        d1_l % Left boundary first derivative
        d1_r % Right boundary first derivative
        m % Number of grid points.
        h % Step size
        x % grid
        borrowing % Struct with borrowing limits for different norm matrices
    end

    methods
        function obj = D2Variable(m,lim,order)

            x_l = lim{1};
            x_r = lim{2};
            L = x_r-x_l;
            obj.h = L/(m-1);
            obj.x = linspace(x_l,x_r,m)';

            switch order
                case 4
                    [obj.H, obj.HI, obj.D1, obj.D2, obj.e_l,...
                        obj.e_r, obj.d1_l, obj.d1_r] = ...
                        sbp.implementations.d2_variable_4(m,obj.h);
                    obj.borrowing.M.d1 = 0.2505765857;
                otherwise
                    error('Invalid operator order %d.',order);
            end

            obj.m = m;
            obj.M = [];

        end

        function str = string(obj)
            str = [class(obj) '_' num2str(obj.order)];
        end
    end


end





