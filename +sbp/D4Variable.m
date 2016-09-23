classdef D4CompatibleVariable < sbp.OpSet
    properties
        D1   % SBP operator approximating first derivative
        H    % Norm matrix
        HI   % H^-1
        Q    % Skew-symmetric matrix
        e_l  % Left boundary operator
        e_r  % Right boundary operator
        D2   % SBP operator for second derivative
        M    % Norm matrix, second derivative
        d1_l % Left boundary first derivative
        d1_r % Right boundary first derivative
        D3   % SBP operator for third derivative
        Q3   % Skew-symmetric matrix in third derivative
        d2_l % Left boundary second derivative
        d2_r % Right boundary second derivative
        D4   % SBP operator for fourth derivative
        M4   % Norm matrix, fourth derivative
        d3_l % Left boundary third derivative
        d3_r % Right boundary third derivative
        m    % Number of grid points.
        h    % Step size
        x    % grid
        borrowing % Struct with borrowing limits for different norm matrices
    end

    methods
        function obj = D4CompatibleVariable(m, lim, order, opt)
            default_arg('opt', '')

            x_l = lim{1};
            x_r = lim{2};
            L = x_r-x_l;
            obj.h = L/(m-1);
            obj.x = linspace(x_l, x_r,m)';

            if order == 2
                [obj.H, obj.HI, ~, obj.D2, ~, obj.D4, obj.e_l, obj.e_r,...
                 obj.M4, ~, obj.d2_l, obj.d2_r, obj.d3_l,...
                    obj.d3_r, obj.d1_l, obj.d1_r] =...
                    sbp.implementations.d4_compatible_halfvariable_2(m,obj.h);
                obj.borrowing.N.S2 = 1.2500;
                obj.borrowing.N.S3 = 0.4000;

            elseif order == 4
                switch opt
                    case 'min_boundary_points'
                        [H, HI, D2, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = ...
                            sbp.higher_variable4_min_boundary_points(m, obj.h);
                        % obj.borrowing.N.S2 = 0.5055;
                        % obj.borrowing.N.S3 = 0.9290;
                    otherwise
                        [H, HI, D2, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = ...
                            sbp.implementations.d4_compatible_halfvariable_4(m, obj.h);
                        obj.borrowing.N.S2 = 0.5055;
                        obj.borrowing.N.S3 = 0.9290;
                end

            elseif order == 6
                switch opt
                    case '2'
                        [H, HI, D2, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = ...
                            sbp.higher_variable6_2(m, obj.h);
                        % obj.borrowing.N.S2 = 0.3259;
                        % obj.borrowing.N.S3 = 0.1580;
                    case '3'
                        [H, HI, D2, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = ...
                            sbp.higher_variable6_3(m, obj.h);
                        % obj.borrowing.N.S2 = 0.3259;
                        % obj.borrowing.N.S3 = 0.1580;
                    case 'min_boundary_points'
                        [H, HI, D2, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = ...
                            sbp.higher_variable6_min_boundary_points(m, obj.h);
                        % obj.borrowing.N.S2 = 0.3259;
                        % obj.borrowing.N.S3 = 0.1580;
                    otherwise
                        [H, HI, D2, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = ...
                            sbp.implementations.d4_compatible_halfvariable_6(m, obj.h);
                        obj.borrowing.N.S2 = 0.3259;
                        obj.borrowing.N.S3 = 0.1580;
                end

            elseif order == 8
                switch opt
                    case 'min_boundary_points'
                        [H, HI, D2, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = ...
                            sbp.higher_variable8_min_boundary_points(m, obj.h);
                        % obj.borrowing.N.S2 = 0.3259;
                        % obj.borrowing.N.S3 = 0.1580;
                    otherwise
                        [H, HI, D2, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = ...
                            sbp.higher_variable8_higher_boundary_order(m, obj.h);
                        % obj.borrowing.N.S2 = 0.3259;
                        % obj.borrowing.N.S3 = 0.1580;
                    end
            else
                error('Invalid operator order.');
            end

            obj.m = m;

            obj.H    = H;
            obj.HI   = HI;
            obj.D2   = D2;
            obj.D4   = D4;
            obj.M4   = M4;
            obj.e_l  = e_l;
            obj.e_r  = e_r;
            obj.d1_l = d1_l;
            obj.d1_r = d1_r ;
            obj.d2_l = d2_l;
            obj.d2_r = d2_r;
            obj.d3_l = d3_l;
            obj.d3_r = d3_r;

            obj.D1 = [];
            obj.D3 = [];
        end
    end
end