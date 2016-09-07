classdef D1_nonequidistant_accurate < sbp.OpSet
    properties
        norms % Struct containing norm matrices such as H,Q, M
        boundary  % Struct contanging vectors for boundry point approximations
        derivatives % Struct containging differentiation operators
        borrowing % Struct with borrowing limits for different norm matrices
        m % Number of grid points.
        h % Step size
        x % grid
    end

    methods
        function obj = D1_nonequidistant_accurate(m,L,order)

            if order == 4
                [D1,H,grid,dx] = D1_4th_4BP_2shifts(m,L);
            elseif order == 6
                [D1,H,grid,dx] = D1_6th_6BP_3shifts(m,L);
            elseif order == 8
                [D1,H,grid,dx] = D1_8th_8BP_4shifts(m,L);
            elseif order == 10
                [D1,H,grid,dx] = D1_10th_10BP_5shifts(m,L);
            elseif order == 12
                [D1,H,grid,dx] = D1_12th_12BP_6shifts(m,L);
            else
                error('Invalid operator order %d.',order);
            end

            Q = H*D1;
            e_1 = sparse(m,1);
            e_m = sparse(m,1);
            e_1(1) = 1;
            e_m(m) = 1;
            
            obj.h = dx;
            obj.m = m;
            obj.x = grid;

            obj.norms.H = H;
            obj.norms.HI = HI;
            obj.norms.Q = Q;

            obj.boundary.e_1 = e_1;
            obj.boundary.e_m = e_m;

            obj.derivatives.D1 = D1;
        end
    end

    methods (Static)
        function lambda = smallestGrid(obj)
            error('Not implmented')
        end
    end
end





