classdef InterpOpsOP < sbp.InterpOps
    properties

        % Structs of interpolation operators, fields .good and .bad
        Iu2v
        Iv2u
    end

    methods
        % m_u, m_v         --   number of grid points along the interface
        % order_u, order_v --   order of accuracy in the different blocks
        function obj = InterpOpsOP(m_u, m_v, order_u, order_v)

            assert(order_u == order_v,...
                    'InterpOpsOP: Different orders of accuracy not available');

            switch order_u
            case 2
                intOpSet = @sbp.implementations.intOpOP_orders_2to2_ratio2to1;
            case 4
                intOpSet = @sbp.implementations.intOpOP_orders_4to4_ratio2to1;
            case 6
                intOpSet = @sbp.implementations.intOpOP_orders_6to6_ratio2to1;
            case 8
                intOpSet = @sbp.implementations.intOpOP_orders_8to8_ratio2to1;
            otherwise
                error('InterpOpsOP: Order of accuracy %d not available.', order_u);
            end

            Iu2v = struct;
            Iv2u = struct;

            if (m_u-1)/(m_v-1) == 2
                % Block u is fine, v is coarse
                m_C = m_v;
                [Iv2u.good, Iu2v.bad] = intOpSet(m_C, 1, 'C2F');
                [Iv2u.bad, Iu2v.good] = intOpSet(m_C, 1, 'F2C');

            elseif (m_v-1)/(m_u-1) == 2
                % Block v is fine, u is coarse
                m_C = m_u;
                [Iu2v.good, Iv2u.bad] = intOpSet(m_C, 1, 'C2F');
                [Iu2v.bad, Iv2u.good] = intOpSet(m_C, 1, 'F2C');
            else
                error('InterpOpsOP: Interpolation operators for grid ratio %f have not yet been constructed', (m_u-1)/(m_v-1));
            end

            obj.Iu2v = Iu2v;
            obj.Iv2u = Iv2u;

        end

        function str = string(obj)
            str = [class(obj)];
        end

    end
end
