classdef InterpAWW < sbp.InterpOps
    properties

        % Interpolation operators
        IC2F
        IF2C

        % Orders used on coarse and fine sides
        order_C
        order_F

        % Grid points, refinement ratio.
        ratio
        m_C
        m_F

        % Boundary accuracy of IC2F and IF2C.
        acc_C2F
        acc_F2C
    end

    methods
        % accOp : String, 'C2F' or 'F2C'. Specifies which of the operators
        % should have higher accuracy.
        function obj = InterpAWW(m_C,m_F,order_C,order_F,accOp)

            ratio = (m_F-1)/(m_C-1);
            h_C = 1;

            assert(order_C == order_F,...
                    'Error: Different orders of accuracy not available');

            switch ratio
            case 2
                switch order_C
                case 2
                    [IC2F,IF2C] = sbp.implementations.intOpAWW_orders_2to2_ratio2to1(m_C, h_C, accOp);
                case 4
                    [IC2F,IF2C] = sbp.implementations.intOpAWW_orders_4to4_ratio2to1(m_C, h_C, accOp);
                case 6
                    [IC2F,IF2C] = sbp.implementations.intOpAWW_orders_6to6_ratio2to1(m_C, h_C, accOp);
                case 8
                    [IC2F,IF2C] = sbp.implementations.intOpAWW_orders_8to8_ratio2to1(m_C, h_C, accOp);
                otherwise
                    error(['Order ' num2str(order_C) ' not available.']);
                end
            otherwise
                error(['Grid ratio ' num2str(ratio) ' not available']);
            end

            obj.IC2F = IC2F;
            obj.IF2C = IF2C;
            obj.order_C = order_C;
            obj.order_F = order_F;
            obj.ratio = ratio;
            obj.m_C = m_C;
            obj.m_F = m_F;

        end

        function str = string(obj)
            str = [class(obj) '_orders' num2str(obj.order_F) 'to'...
            num2str(obj.order_C) '_ratio' num2str(obj.ratio) 'to1'];
        end

    end
end
