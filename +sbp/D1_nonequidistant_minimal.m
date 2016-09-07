classdef D1_nonequidistant_minimal < sbp.OpSet
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
        function obj = D1_nonequidistant_minimal(m,L,order)

            if order == 4
                [D1,H,grid,dx] = D1_minimal_4th_3BP_1shifts(m,L);
            elseif order == 6
                [D1,H,grid,dx] = D1_minimal_6th_5BP_2shifts(m,L);
            elseif order == 8
                [D1,H,grid,dx] = D1_minimal_8th_6BP_2shifts(m,L);
            elseif order == 10
                [D1,H,grid,dx] = D1_minimal_10th_8BP_3shifts(m,L);
            elseif order == 12
                [D1,H,grid,dx] = D1_minimal_12th_10BP_4shifts(m,L);
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





