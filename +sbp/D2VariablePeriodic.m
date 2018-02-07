classdef D2VariablePeriodic < sbp.OpSet
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
        function obj = D2VariablePeriodic(m,lim,order)

            x_l = lim{1};
            x_r = lim{2};
            L = x_r-x_l;
            obj.h = L/m;
            x = linspace(x_l,x_r,m+1)';
            obj.x = x(1:end-1);

            switch order

                case 6
                     error('Not impl')

                    [obj.H, obj.HI, obj.D1, obj.D2, ...
                    ~, obj.e_l, obj.e_r, ~, ~, ~, ~, ~,...
                     obj.d1_l, obj.d1_r] = ...
                        sbp.implementations.d4_variable_periodic_6(m, obj.h);
                    obj.borrowing.M.d1 = 0.1878;
                    obj.borrowing.R.delta_D = 0.3696;
                    % Borrowing e^T*D1 - d1 from R

                case 4
                    error('Not impl')

                    [obj.H, obj.HI, obj.D1, obj.D2, obj.e_l,...
                        obj.e_r, obj.d1_l, obj.d1_r] = ...
                        sbp.implementations.d2_variable_periodic_4(m,obj.h);
                    obj.borrowing.M.d1 = 0.2505765857;

                    obj.borrowing.R.delta_D = 0.577587500088313;
                    % Borrowing e^T*D1 - d1 from R
                case 2
                    [obj.H, obj.HI, obj.D1, obj.D2, obj.e_l,...
                        obj.e_r, obj.d1_l, obj.d1_r] = ...
                        sbp.implementations.d2_variable_periodic_2(m,obj.h);
                    obj.borrowing.M.d1 = 0.3636363636; 
                    % Borrowing const taken from Virta 2014

                    obj.borrowing.R.delta_D = 1.000000538455350;
                    % Borrowing e^T*D1 - d1 from R
                    
                otherwise
                    error('Invalid operator order %d.',order);
            end
            obj.borrowing.H11 = obj.H(1,1)/obj.h; % First element in H/h,
            obj.m = m;
            obj.M = [];
        end
        function str = string(obj)
            str = [class(obj) '_' num2str(obj.order)];
        end
    end


end
