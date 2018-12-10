% Order-preserving (OP) interpolation operators, see
% Almquist, Wang, Werpers,
% "Order-Preserving Interpolation for Summation-by-Parts Operators
% at Non-Conforming Interfaces", https://arxiv.org/abs/1806.01931
%
% Let ^* denote the adjoint. These operators satsify
%
% Iuv2.good = Iv2u.bad^*
% Iv2u.good = Iu2v.bad^*
%
% The .bad operators have the same order of accuracy as the operators
% by Mattsson and Carpenter (MC) in InterpOpsMC, i.e. order p,
% if the interior stencil is order 2p. The .good operators are
% one order more accurate, i.e. order p+1.
%
% For PDEs of second order in space, the OP operators allow for the same
% convergence rate as with conforming interfaces, which is an improvement
% by one order compared what is possible with the MC operators.
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
