classdef D1Gauss < sbp.OpSet
    % Diagonal-norm SBP operators based on the Gauss quadrature formula
    % with m nodes, which is of degree 2m-1. Hence, The operator D1 is
    % accurate of order m.
    properties
        D1 % SBP operator approximating first derivative
        H % Norm matrix
        HI % H^-1
        Q % Skew-symmetric matrix
        e_l % Left boundary operator
        e_r % Right boundary operator
        m % Number of grid points.
        h % Step size
        x % grid
        borrowing % Struct with borrowing limits for different norm matrices
    end

    methods
        function obj = D1Gauss(m,lim)

            x_l = lim{1};
            x_r = lim{2};
            L = x_r-x_l;

            switch m
                case 4
                    [obj.D1,obj.H,obj.x,obj.h,obj.e_l,obj.e_r] = ...
                        sbp.implementations.d1_gauss_4(L);
                otherwise
                    error('Invalid number of points: %d.', m);
            end


            obj.x = obj.x + x_l;
            obj.HI = inv(obj.H);
            obj.Q = obj.H*obj.D1 - obj.e_r*obj.e_r' + obj.e_l*obj.e_l';

            obj.borrowing = [];
        end

        function str = string(obj)
            str = [class(obj) '_' num2str(obj.order)];
        end
    end
end
