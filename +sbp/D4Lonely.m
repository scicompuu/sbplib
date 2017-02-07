classdef D4Lonely < sbp.OpSet
    properties
        m    % Number of grid points.
        h    % Step size
        x    % grid
        H    % Norm matrix
        HI   % H^-1
        D4   % SBP operator for fourth derivative
        M4   % Norm matrix, fourth derivative
        e_l,  e_r  % Left and right boundary operator
        d1_l, d1_r % Left and right boundary first derivative
        d2_l, d2_r % Left and right boundary second derivative
        d3_l, d3_r % Left and right boundary third derivative
        borrowing % Struct with borrowing limits for different norm matrices
        opt
        order
    end

    methods
        function obj = D4Lonely(m, lim, order, opt)
            default_arg('opt', '')

            obj.opt = opt;
            obj.order = order;

            x_l = lim{1};
            x_r = lim{2};
            L = x_r-x_l;
            obj.h = L/(m-1);
            obj.x = linspace(x_l, x_r,m)';

            if order == 2
                switch opt
                    case ''
                        [H, HI, D1, D2, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = ...
                            sbp.implementations.d4_variable_2(m, obj.h);
                        obj.borrowing.N.S2 = 1.2500;
                        obj.borrowing.N.S3 = 0.4000;
                    otherwise
                        error('Invalid operator option.');
                end

            elseif order == 4
                switch opt
                    case 'min_boundary_points'
                        [H, HI, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = ...
                            sbp.implementations.d4_lonely_4_min_boundary_points(m, obj.h);
                        obj.borrowing.N.S2 = 0.6244;
                        obj.borrowing.N.S3 = 1.3961;
                    case ''
                        [H, HI, D1, D2, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = ...
                            sbp.implementations.d4_variable_4(m, obj.h);
                        obj.borrowing.N.S2 = 0.5055;
                        obj.borrowing.N.S3 = 0.9290;
                    otherwise
                        error('Invalid operator option.');
                end

            elseif order == 6
                switch opt
                    case '2'
                        [H, HI, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = ...
                            sbp.implementations.d4_lonely_6_2(m, obj.h);
                        obj.borrowing.N.S2 = 0.2931;
                        obj.borrowing.N.S3 = 0.0807;
                    case '3'
                        [H, HI, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = ...
                            sbp.implementations.d4_lonely_6_3(m, obj.h);
                        obj.borrowing.N.S2 = 0.2842;
                        obj.borrowing.N.S3 = 0.0709;
                    case 'min_boundary_points'
                        [H, HI, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = ...
                            sbp.implementations.d4_lonely_6_min_boundary_points(m, obj.h);
                        obj.borrowing.N.S2 = 0.3569;
                        obj.borrowing.N.S3 = 0.1908;
                    case ''
                        [H, HI, D1, D2, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = ...
                            sbp.implementations.d4_variable_6(m, obj.h);
                        obj.borrowing.N.S2 = 0.3259;
                        obj.borrowing.N.S3 = 0.1580;
                    otherwise
                        error('Invalid operator option.');
                end

            elseif order == 8
                switch opt
                    case 'min_boundary_points'
                        [H, HI, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = ...
                            sbp.implementations.d4_lonely_8_min_boundary_points(m, obj.h);
                        obj.borrowing.N.S2 = 0.2804;
                        obj.borrowing.N.S3 = 0.0740;
                    case ''
                        [H, HI, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = ...
                            sbp.implementations.d4_lonely_8_higher_boundary_order(m, obj.h);
                        obj.borrowing.N.S2 = 0.2475;
                        obj.borrowing.N.S3 = 0.0401;
                    otherwise
                        error('Invalid operator option.');
                    end
            else
                error('Invalid operator order.');
            end

            obj.m = m;

            obj.H    = H;
            obj.HI   = HI;
            obj.D4   = D4;
            obj.M4   = M4;
            obj.e_l  = e_l;
            obj.e_r  = e_r;
            obj.d1_l = d1_l;
            obj.d1_r = d1_r;
            obj.d2_l = d2_l;
            obj.d2_r = d2_r;
            obj.d3_l = d3_l;
            obj.d3_r = d3_r;
        end

        function str = string(obj)
            str = [class(obj) '_' num2str(obj.order) '_' obj.opt];
        end
    end
end