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
        borrowing % Struct with borrowing limits for different norm matrices
    end
    
    methods
        function obj = D1Nonequidistant(m,L,order,option)
            
            default_arg('option','Accurate');
            % 'Accurate' operators are optimized for accuracy
            % 'Minimal' operators have the smallest possible boundary
            %  closure
            
            switch option
                
                case {'Accurate','accurate','A'}
                    
                    if order == 4
                        [obj.D1,obj.H,obj.x,obj.h] = ...
                            sbp.implementations.d1_noneq_4(m,L);
                    elseif order == 6
                        [obj.D1,obj.H,obj.x,obj.h] = ...
                            sbp.implementations.d1_noneq_6(m,L);
                    elseif order == 8
                        [obj.D1,obj.H,obj.x,obj.h] = ...
                            sbp.implementations.d1_noneq_8(m,L);
                    elseif order == 10
                        [obj.D1,obj.H,obj.x,obj.h] = ...
                            sbp.implementations.d1_noneq_10(m,L);
                    elseif order == 12
                        [obj.D1,obj.H,obj.x,obj.h] = ...
                            sbp.implementations.d1_noneq_12(m,L);
                    else
                        error('Invalid operator order %d.',order);
                    end
                    
                case {'Minimal','minimal','M'}
                    
                    if order == 4
                        [obj.D1,obj.H,obj.x,obj.h] = ...
                            sbp.implementations.d1_noneq_minimal_4(m,L);
                    elseif order == 6
                        [obj.D1,obj.H,obj.x,obj.h] = ...
                            sbp.implementations.d1_noneq_minimal_6(m,L);
                    elseif order == 8
                        [obj.D1,obj.H,obj.x,obj.h] = ...
                            sbp.implementations.d1_noneq_minimal_8(m,L);
                    elseif order == 10
                        [obj.D1,obj.H,obj.x,obj.h] = ...
                            sbp.implementations.d1_noneq_minimal_10(m,L);
                    elseif order == 12
                        [obj.D1,obj.H,obj.x,obj.h] = ...
                            sbp.implementations.d1_noneq_minimal_12(m,L);
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
            
            obj.borrowing = [];
            
        end
    end
    
    methods (Static)
        function lambda = smallestGrid(obj)
            error('Not implmented')
        end
    end
end





