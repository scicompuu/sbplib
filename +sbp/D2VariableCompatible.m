classdef D2VariableCompatible < sbp.OpSet
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
        function obj = D2VariableCompatible(m,lim,order)

            x_l = lim{1};
            x_r = lim{2};
            L = x_r-x_l;
            obj.h = L/(m-1);
            obj.x = linspace(x_l,x_r,m)';

            switch order

                case 6

                    [obj.H, obj.HI, obj.D1, D2, ...
                    ~, obj.e_l, obj.e_r, ~, ~, ~, ~, ~,...
                     d1_l, d1_r] = ...
                        sbp.implementations.d4_variable_6(m, obj.h);

                case 4
                    [obj.H, obj.HI, obj.D1, D2, obj.e_l,...
                        obj.e_r, d1_l, d1_r] = ...
                        sbp.implementations.d2_variable_4(m,obj.h);
                case 2
                    [obj.H, obj.HI, obj.D1, D2, obj.e_l,...
                        obj.e_r, d1_l, d1_r] = ...
                        sbp.implementations.d2_variable_2(m,obj.h);

                otherwise
                    error('Invalid operator order %d.',order);
            end
            obj.borrowing.H11 = obj.H(1,1)/obj.h; % First element in H/h,
            obj.borrowing.M.d1 = obj.H(1,1)/obj.h; % First element in H/h is borrowing also for M
            obj.borrowing.R.delta_D = inf;
            obj.m = m;
            obj.M = [];


            D1 = obj.D1;
            e_r = obj.e_r;
            e_l = obj.e_l;

            % D2 = Hinv * (-M + br*er*d1r^T - bl*el*d1l^T);
            % Replace d1' by e'*D1 in D2.
            correction_l = obj.HI*(e_l*d1_l' - e_l*e_l'*D1);
            correction_r = - obj.HI*(e_r*d1_r' - e_r*e_r'*D1);

            D2_compatible = @(b) D2(b) + b(1)*correction_l + b(m)*correction_r;

            obj.D2 = D2_compatible;
            obj.d1_l = (e_l'*D1)';
            obj.d1_r = (e_r'*D1)';

        end
        function str = string(obj)
            str = [class(obj) '_' num2str(obj.order)];
        end
    end


end





