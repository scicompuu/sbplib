classdef D1Upwind < sbp.OpSet
    properties
        D1 % SBP operator approximating first derivative
        H % Norm matrix
        HI % H^-1
        Q % Skew-symmetric matrix
        e_1 % Left boundary operator
        e_m % Right boundary operator
        m % Number of grid points.
        h % Step size
        x % grid
        borrowing % Struct with borrowing limits for different norm matrices
    end

    methods
        function obj = D1Upwind(m,lim,order)
            
            x_l = lim{1};
            x_r = lim{2};
            L = x_r-x_l;
            obj.h = L/(m-1);
            obj.x = linspace(x_l,x_r,m)';

            switch order
                case 2
                    [obj.H, obj.HI, obj.Dp, obj.Dm, obj.e_1, obj.e_m] = ...
                        sbp.implementations.d1_upwind_2(m,obj.h);
                case 3
                    [obj.H, obj.HI, obj.Dp, obj.Dm, obj.e_1, obj.e_m] = ...
                        sbp.implementations.d1_upwind_3(m,obj.h);
                case 4
                    [obj.H, obj.HI, obj.Dp, obj.Dm, obj.e_1, obj.e_m] = ...
                        sbp.implementations.d1_upwind_4(m,obj.h);
                case 5
                    [obj.H, obj.HI, obj.Dp, obj.Dm, obj.e_1, obj.e_m] = ...
                        sbp.implementations.d1_upwind_5(m,obj.h);
                case 6
                    [obj.H, obj.HI, obj.Dp, obj.Dm, obj.e_1, obj.e_m] = ...
                        sbp.implementations.d1_upwind_6(m,obj.h);
                case 7
                    [obj.H, obj.HI, obj.Dp, obj.Dm, obj.e_1, obj.e_m] = ...
                        sbp.implementations.d1_upwind_7(m,obj.h);
                case 8
                    [obj.H, obj.HI, obj.Dp, obj.Dm, obj.e_1, obj.e_m] = ...
                        sbp.implementations.d1_upwind_8(m,obj.h);
                case 9
                    [obj.H, obj.HI, obj.Dp, obj.Dm, obj.e_1, obj.e_m] = ...
                        sbp.implementations.d1_upwind_9(m,obj.h);
                otherwise
                    error('Invalid operator order %d.',order);
            end

            obj.m = m;
        	obj.borrowing = [];

        end
    end


end





