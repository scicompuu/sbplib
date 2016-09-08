classdef D1Nonequidistant < sbp.OpSet
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
    end
    
    methods
        function obj = D1Nonequidistant(m,L,order,option)
            
            default_arg('option','Accurate');
            % 'Accurate' operators are optimized for accuracy
            % 'Minimal' operators have the smallest possible boundary
            % closure
            
            switch option
                
                case {'Accurate','accurate','A'}
                    
                    if order == 4
                        [obj.D1,obj.H,obj.x,obj.h] = sbp.D1_4th_4BP_2shifts(m,L);
                    elseif order == 6
                        [obj.D1,obj.H,obj.x,obj.h] = sbp.D1_6th_6BP_3shifts(m,L);
                    elseif order == 8
                        [obj.D1,obj.H,obj.x,obj.h] = sbp.D1_8th_8BP_4shifts(m,L);
                    elseif order == 10
                        [obj.D1,obj.H,obj.x,obj.h] = sbp.D1_10th_10BP_5shifts(m,L);
                    elseif order == 12
                        [obj.D1,obj.H,obj.x,obj.h] = sbp.D1_12th_12BP_6shifts(m,L);
                    else
                        error('Invalid operator order %d.',order);
                    end
                    
                case {'Minimal','minimal','M'}
                    
                    if order == 4
                        [obj.D1,obj.H,obj.x,obj.h] = sbp.D1_minimal_4th_3BP_1shifts(m,L);
                    elseif order == 6
                        [obj.D1,obj.H,obj.x,obj.h] = sbp.D1_minimal_6th_5BP_2shifts(m,L);
                    elseif order == 8
                        [obj.D1,obj.H,obj.x,obj.h] = sbp.D1_minimal_8th_6BP_2shifts(m,L);
                    elseif order == 10
                        [obj.D1,obj.H,obj.x,obj.h] = sbp.D1_minimal_10th_8BP_3shifts(m,L);
                    elseif order == 12
                        [obj.D1,obj.H,obj.x,obj.h] = sbp.D1_minimal_12th_10BP_4shifts(m,L);
                    else
                        error('Invalid operator order %d.',order);
                    end
                    
            end
            
            obj.e_1 = sparse(m,1);
            obj.e_m = sparse(m,1);
            obj.e_1(1) = 1;
            obj.e_m(m) = 1;
            
            obj.HI = inv(obj.H);
            obj.Q = obj.H*obj.D1 - obj.e_m*obj.e_m' + obj.e_0*obj.e_0';
            
        end
    end
    
    methods (Static)
        function lambda = smallestGrid(obj)
            error('Not implmented')
        end
    end
end





