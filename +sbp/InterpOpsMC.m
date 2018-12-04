classdef InterpOpsMC < sbp.InterpOps
    properties

        % Structs of interpolation operators, fields .good and .bad
        % Here .good and .bad are the same, but this makes them fit in the
        % OP (order-preserving) framework.
        Iu2v
        Iv2u
    end

    methods
        % m_u, m_v         --   number of grid points along the interface
        % order_u, order_v --   order of accuracy in the different blocks
        function obj = InterpOpsMC(m_u, m_v, order_u, order_v)

            assert(order_u == order_v,...
                    'InterpOpsMC: Different orders of accuracy not available');

            switch order_u
            case 2
                intOpSet = @sbp.implementations.intOpMC_orders_2to2_ratio2to1;
            case 4
                intOpSet = @sbp.implementations.intOpMC_orders_4to4_ratio2to1;
            case 6
                intOpSet = @sbp.implementations.intOpMC_orders_6to6_ratio2to1;
            case 8
                intOpSet = @sbp.implementations.intOpMC_orders_8to8_ratio2to1;
            otherwise
                error('InterpOpsMC: Order of accuracy %d not available.', order_u);
            end

            Iu2v = struct;
            Iv2u = struct;

            if (m_u-1)/(m_v-1) == 2
                % Block u is fine, v is coarse
                m_C = m_v;
                [Iv2u.good, Iu2v.bad] = intOpSet(m_C);
                Iv2u.bad = Iv2u.good;
                Iu2v.good = Iu2v.bad;

            elseif (m_v-1)/(m_u-1) == 2
                % Block v is fine, u is coarse
                m_C = m_u;
                [Iu2v.good, Iv2u.bad] = intOpSet(m_C);
                Iu2v.bad = Iu2v.good;
                Iv2u.good = Iv2u.bad;
            else
                error('InterpOpsMC: Interpolation operators for grid ratio %f have not yet been constructed', (m_u-1)/(m_v-1));
            end

            obj.Iu2v = Iu2v;
            obj.Iv2u = Iv2u;

        end

        function str = string(obj)
            str = [class(obj)];
        end

    end
end
