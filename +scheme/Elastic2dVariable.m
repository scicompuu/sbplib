classdef Elastic2dVariable < scheme.Scheme

% Discretizes the elastic wave equation:
% rho u_{i,tt} = di lambda dj u_j + dj mu di u_j + dj mu dj u_i
% opSet should be cell array of opSets, one per dimension. This
% is useful if we have periodic BC in one direction.

    properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing

        grid
        dim

        order % Order of accuracy for the approximation

        % Diagonal matrices for varible coefficients
        LAMBDA % Variable coefficient, related to dilation
        MU     % Shear modulus, variable coefficient
        RHO, RHOi % Density, variable

        D % Total operator
        D1 % First derivatives

        % Second derivatives
        D2_lambda
        D2_mu

        % Traction operators used for BC
        T_l, T_r
        tau_l, tau_r

        H, Hi, H_1D % Inner products
        e_l, e_r


        d1_l, d1_r % Normal derivatives at the boundary
        E % E{i}^T picks out component i

        H_boundary % Boundary inner products

        % Kroneckered norms and coefficients
        RHOi_kron
        Hi_kron

        % Borrowing constants of the form gamma*h, where gamma is a dimensionless constant.
        theta_R % Borrowing (d1- D1)^2 from R
        theta_H % First entry in norm matrix
        theta_M % Borrowing d1^2 from M.

        % Structures used for adjoint optimization
        B
    end

    methods

        % The coefficients can either be function handles or grid functions
        function obj = Elastic2dVariable(g ,order, lambda, mu, rho, opSet)
            default_arg('opSet',{@sbp.D2Variable, @sbp.D2Variable});
            default_arg('lambda', @(x,y) 0*x+1);
            default_arg('mu', @(x,y) 0*x+1);
            default_arg('rho', @(x,y) 0*x+1);
            dim = 2;

            assert(isa(g, 'grid.Cartesian'))

            if isa(lambda, 'function_handle')
                lambda = grid.evalOn(g, lambda);
            end
            if isa(mu, 'function_handle')
                mu = grid.evalOn(g, mu);
            end
            if isa(rho, 'function_handle')
                rho = grid.evalOn(g, rho);
            end

            m = g.size();
            m_tot = g.N();

            h = g.scaling();
            lim = g.lim;
            if isempty(lim)
                x = g.x;
                lim = cell(length(x),1);
                for i = 1:length(x)
                    lim{i} = {min(x{i}), max(x{i})};
                end
            end

            % 1D operators
            ops = cell(dim,1);
            for i = 1:dim
                ops{i} = opSet{i}(m(i), lim{i}, order);
            end

            % Borrowing constants
            for i = 1:dim
                obj.theta_R{i} = h(i)*ops{i}.borrowing.R.delta_D;
                obj.theta_H{i} = h(i)*ops{i}.borrowing.H11;
                obj.theta_M{i} = h(i)*ops{i}.borrowing.M.d1;
            end

            I = cell(dim,1);
            D1 = cell(dim,1);
            D2 = cell(dim,1);
            H = cell(dim,1);
            Hi = cell(dim,1);
            e_l = cell(dim,1);
            e_r = cell(dim,1);
            d1_l = cell(dim,1);
            d1_r = cell(dim,1);

            for i = 1:dim
                I{i} = speye(m(i));
                D1{i} = ops{i}.D1;
                D2{i} = ops{i}.D2;
                H{i} =  ops{i}.H;
                Hi{i} = ops{i}.HI;
                e_l{i} = ops{i}.e_l;
                e_r{i} = ops{i}.e_r;
                d1_l{i} = ops{i}.d1_l;
                d1_r{i} = ops{i}.d1_r;
            end

            %====== Assemble full operators ========
            LAMBDA = spdiag(lambda);
            obj.LAMBDA = LAMBDA;
            MU = spdiag(mu);
            obj.MU = MU;
            RHO = spdiag(rho);
            obj.RHO = RHO;
            obj.RHOi = inv(RHO);

            obj.D1 = cell(dim,1);
            obj.D2_lambda = cell(dim,1);
            obj.D2_mu = cell(dim,1);
            obj.e_l = cell(dim,1);
            obj.e_r = cell(dim,1);
            obj.d1_l = cell(dim,1);
            obj.d1_r = cell(dim,1);

            % D1
            obj.D1{1} = kron(D1{1},I{2});
            obj.D1{2} = kron(I{1},D1{2});

            % Boundary operators
            obj.e_l{1} = kron(e_l{1},I{2});
            obj.e_l{2} = kron(I{1},e_l{2});
            obj.e_r{1} = kron(e_r{1},I{2});
            obj.e_r{2} = kron(I{1},e_r{2});

            obj.d1_l{1} = kron(d1_l{1},I{2});
            obj.d1_l{2} = kron(I{1},d1_l{2});
            obj.d1_r{1} = kron(d1_r{1},I{2});
            obj.d1_r{2} = kron(I{1},d1_r{2});

            % D2
            for i = 1:dim
                obj.D2_lambda{i} = sparse(m_tot);
                obj.D2_mu{i} = sparse(m_tot);
            end
            ind = grid.funcToMatrix(g, 1:m_tot);

            for i = 1:m(2)
                D_lambda = D2{1}(lambda(ind(:,i)));
                D_mu = D2{1}(mu(ind(:,i)));

                p = ind(:,i);
                obj.D2_lambda{1}(p,p) = D_lambda;
                obj.D2_mu{1}(p,p) = D_mu;
            end

            for i = 1:m(1)
                D_lambda = D2{2}(lambda(ind(i,:)));
                D_mu = D2{2}(mu(ind(i,:)));

                p = ind(i,:);
                obj.D2_lambda{2}(p,p) = D_lambda;
                obj.D2_mu{2}(p,p) = D_mu;
            end

            % Quadratures
            obj.H = kron(H{1},H{2});
            obj.Hi = inv(obj.H);
            obj.H_boundary = cell(dim,1);
            obj.H_boundary{1} = H{2};
            obj.H_boundary{2} = H{1};
            obj.H_1D = {H{1}, H{2}};

            % E{i}^T picks out component i.
            E = cell(dim,1);
            I = speye(m_tot,m_tot);
            for i = 1:dim
                e = sparse(dim,1);
                e(i) = 1;
                E{i} = kron(I,e);
            end
            obj.E = E;

            % Differentiation matrix D (without SAT)
            D2_lambda = obj.D2_lambda;
            D2_mu = obj.D2_mu;
            D1 = obj.D1;
            D = sparse(dim*m_tot,dim*m_tot);
            d = @kroneckerDelta;    % Kronecker delta
            db = @(i,j) 1-d(i,j); % Logical not of Kronecker delta
            for i = 1:dim
                for j = 1:dim
                    D = D + E{i}*inv(RHO)*( d(i,j)*D2_lambda{i}*E{j}' +...
                                            db(i,j)*D1{i}*LAMBDA*D1{j}*E{j}' ...
                                          );
                    D = D + E{i}*inv(RHO)*( d(i,j)*D2_mu{i}*E{j}' +...
                                            db(i,j)*D1{j}*MU*D1{i}*E{j}' + ...
                                            D2_mu{j}*E{i}' ...
                                          );
                end
            end
            obj.D = D;
            %=========================================%'

            % Numerical traction operators for BC.
            % Because d1 =/= e0^T*D1, the numerical tractions are different
            % at every boundary.
            T_l = cell(dim,1);
            T_r = cell(dim,1);
            tau_l = cell(dim,1);
            tau_r = cell(dim,1);
            % tau^{j}_i = sum_k T^{j}_{ik} u_k

            d1_l = obj.d1_l;
            d1_r = obj.d1_r;
            e_l = obj.e_l;
            e_r = obj.e_r;
            D1 = obj.D1;

            % Loop over boundaries
            for j = 1:dim
                T_l{j} = cell(dim,dim);
                T_r{j} = cell(dim,dim);
                tau_l{j} = cell(dim,1);
                tau_r{j} = cell(dim,1);

                LAMBDA_l = e_l{j}'*LAMBDA*e_l{j};
                LAMBDA_r = e_r{j}'*LAMBDA*e_r{j};
                MU_l = e_l{j}'*MU*e_l{j};
                MU_r = e_r{j}'*MU*e_r{j};

                [~, n_l] = size(e_l{j});
                [~, n_r] = size(e_r{j});

                % Loop over components
                for i = 1:dim
                    tau_l{j}{i} = sparse(n_l, dim*m_tot);
                    tau_r{j}{i} = sparse(n_r, dim*m_tot);
                    for k = 1:dim
                        T_l{j}{i,k} = ...
                        -d(i,j)*LAMBDA_l*(d(i,k)*d1_l{j}' + db(i,k)*e_l{j}'*D1{k})...
                        -d(j,k)*MU_l*(d(i,j)*d1_l{j}' + db(i,j)*e_l{j}'*D1{i})...
                        -d(i,k)*MU_l*d1_l{j}';

                        T_r{j}{i,k} = ...
                        d(i,j)*LAMBDA_r*(d(i,k)*d1_r{j}' + db(i,k)*e_r{j}'*D1{k})...
                        +d(j,k)*MU_r*(d(i,j)*d1_r{j}' + db(i,j)*e_r{j}'*D1{i})...
                        +d(i,k)*MU_r*d1_r{j}';

                        tau_l{j}{i} = tau_l{j}{i} + T_l{j}{i,k}*E{k}';
                        tau_r{j}{i} = tau_r{j}{i} + T_r{j}{i,k}*E{k}';
                    end

                end
            end

            % Transpose T and tau to match boundary operator convention
            for i = 1:dim
                for j = 1:dim
                    tau_l{i}{j} = transpose(tau_l{i}{j});
                    tau_r{i}{j} = transpose(tau_r{i}{j});
                    for k = 1:dim
                        T_l{i}{j,k} = transpose(T_l{i}{j,k});
                        T_r{i}{j,k} = transpose(T_r{i}{j,k});
                    end
                end
            end

            obj.T_l = T_l;
            obj.T_r = T_r;
            obj.tau_l = tau_l;
            obj.tau_r = tau_r;

            % Kroneckered norms and coefficients
            I_dim = speye(dim);
            obj.RHOi_kron = kron(obj.RHOi, I_dim);
            obj.Hi_kron = kron(obj.Hi, I_dim);

            % Misc.
            obj.m = m;
            obj.h = h;
            obj.order = order;
            obj.grid = g;
            obj.dim = dim;

            % B, used for adjoint optimization
            B = cell(dim, 1);
            for i = 1:dim
                B{i} = cell(m_tot, 1);
            end

            for i = 1:dim
                for j = 1:m_tot
                    B{i}{j} = sparse(m_tot, m_tot);
                end
            end

            ind = grid.funcToMatrix(g, 1:m_tot);

            % Direction 1
            for k = 1:m(1)
                c = sparse(m(1),1);
                c(k) = 1;
                [~, B_1D] = ops{1}.D2(c);
                for l = 1:m(2)
                    p = ind(:,l);
                    B{1}{(k-1)*m(2) + l}(p, p) = B_1D;
                end
            end

            % Direction 2
            for k = 1:m(2)
                c = sparse(m(2),1);
                c(k) = 1;
                [~, B_1D] = ops{2}.D2(c);
                for l = 1:m(1)
                    p = ind(l,:);
                    B{2}{(l-1)*m(2) + k}(p, p) = B_1D;
                end
            end

            obj.B = B;

        end


        % Closure functions return the operators applied to the own domain to close the boundary
        % Penalty functions return the operators to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       bc                  is a cell array of component and bc type, e.g. {1, 'd'} for Dirichlet condition
        %                           on the first component.
        %       data                is a function returning the data that should be applied at the boundary.
        %       neighbour_scheme    is an instance of Scheme that should be interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to interface to.
        function [closure, penalty] = boundary_condition(obj, boundary, bc, tuning)
            default_arg('tuning', 1.2);

            assert( iscell(bc), 'The BC type must be a 2x1 cell array' );
            comp = bc{1};
            type = bc{2};

            % j is the coordinate direction of the boundary
            j = obj.get_boundary_number(boundary);
            [e, T, tau, H_gamma] = obj.getBoundaryOperator({'e','T','tau','H'}, boundary);


            E = obj.E;
            Hi = obj.Hi;
            LAMBDA = obj.LAMBDA;
            MU = obj.MU;
            RHOi = obj.RHOi;

            dim = obj.dim;
            m_tot = obj.grid.N();

            % Preallocate
            closure = sparse(dim*m_tot, dim*m_tot);
            penalty = sparse(dim*m_tot, m_tot/obj.m(j));

            k = comp;
            switch type

            % Dirichlet boundary condition
            case {'D','d','dirichlet','Dirichlet'}

                alpha = obj.getBoundaryOperator('alpha', boundary);

                % Loop over components that Dirichlet penalties end up on
                for i = 1:dim
                    C = transpose(T{k,i});
                    A = -tuning*e*transpose(alpha{i,k});
                    B = A + e*C;
                    closure = closure + E{i}*RHOi*Hi*B'*e*H_gamma*(e'*E{k}' );
                    penalty = penalty - E{i}*RHOi*Hi*B'*e*H_gamma;
                end

            % Free boundary condition
            case {'F','f','Free','free','traction','Traction','t','T'}
                    closure = closure - E{k}*RHOi*Hi*e*H_gamma*tau{k}';
                    penalty = penalty + E{k}*RHOi*Hi*e*H_gamma;

            % Unknown boundary condition
            otherwise
                error('No such boundary condition: type = %s',type);
            end
        end

        % type     Struct that specifies the interface coupling.
        %          Fields:
        %          -- tuning:           penalty strength, defaults to 1.2
        %          -- interpolation:    type of interpolation, default 'none'
        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary,type)

            defaultType.tuning = 1.2;
            defaultType.interpolation = 'none';
            default_struct('type', defaultType);

            switch type.interpolation
            case {'none', ''}
                [closure, penalty] = interfaceStandard(obj,boundary,neighbour_scheme,neighbour_boundary,type);
            case {'op','OP'}
                [closure, penalty] = interfaceNonConforming(obj,boundary,neighbour_scheme,neighbour_boundary,type);
            otherwise
                error('Unknown type of interpolation: %s ', type.interpolation);
            end
        end

        function [closure, penalty] = interfaceStandard(obj,boundary,neighbour_scheme,neighbour_boundary,type)
            tuning = type.tuning;

            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            % Operators without subscripts are from the own domain.

            % Get boundary operators
            e = obj.getBoundaryOperator('e_tot', boundary);
            tau = obj.getBoundaryOperator('tau_tot', boundary);

            e_v = neighbour_scheme.getBoundaryOperator('e_tot', neighbour_boundary);
            tau_v = neighbour_scheme.getBoundaryOperator('tau_tot', neighbour_boundary);

            H_gamma = obj.getBoundaryQuadrature(boundary);

            % Operators and quantities that correspond to the own domain only
            Hi = obj.Hi_kron;
            RHOi = obj.RHOi_kron;

            % Penalty strength operators
            alpha_u = 1/4*tuning*obj.getBoundaryOperator('alpha_tot', boundary);
            alpha_v = 1/4*tuning*neighbour_scheme.getBoundaryOperator('alpha_tot', neighbour_boundary);

            closure = -RHOi*Hi*e*H_gamma*(alpha_u' + alpha_v'*e_v*e');
            penalty = RHOi*Hi*e*H_gamma*(alpha_u'*e*e_v' + alpha_v');

            closure = closure - 1/2*RHOi*Hi*e*H_gamma*tau';
            penalty = penalty - 1/2*RHOi*Hi*e*H_gamma*tau_v';

            closure = closure + 1/2*RHOi*Hi*tau*H_gamma*e';
            penalty = penalty - 1/2*RHOi*Hi*tau*H_gamma*e_v';

        end

        function [closure, penalty] = interfaceNonConforming(obj,boundary,neighbour_scheme,neighbour_boundary,type)
            error('Non-conforming interfaces not implemented yet.');
        end

        % Returns the coordinate number and outward normal component for the boundary specified by the string boundary.
        function [j, nj] = get_boundary_number(obj, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            switch boundary
                case {'w', 'e'}
                    j = 1;
                case {'s', 'n'}
                    j = 2;
            end

            switch boundary
                case {'w', 's'}
                    nj = -1;
                case {'e', 'n'}
                    nj = 1;
            end
        end

        % Returns the boundary operator op for the boundary specified by the string boundary.
        % op -- string
        % Only operators with name *_tot can be used with multiblock.DiffOp.getBoundaryOperator()
        function [varargout] = getBoundaryOperator(obj, op, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})
            assertIsMember(op, {'e', 'e_tot', 'd', 'T', 'tau', 'tau_tot', 'H', 'alpha', 'alpha_tot'})

            switch boundary
                case {'w', 'e'}
                    j = 1;
                case {'s', 'n'}
                    j = 2;
            end

            switch op
                case 'e'
                    switch boundary
                        case {'w', 's'}
                            o = obj.e_l{j};
                        case {'e', 'n'}
                            o = obj.e_r{j};
                    end

                case 'e_tot'
                    e = obj.getBoundaryOperator('e', boundary);
                    I_dim = speye(obj.dim, obj.dim);
                    o = kron(e, I_dim);

                case 'd'
                    switch boundary
                        case {'w', 's'}
                            o = obj.d1_l{j};
                        case {'e', 'n'}
                            o = obj.d1_r{j};
                    end

                case 'T'
                    switch boundary
                        case {'w', 's'}
                            o = obj.T_l{j};
                        case {'e', 'n'}
                            o = obj.T_r{j};
                    end

                case 'tau'
                    switch boundary
                        case {'w', 's'}
                            o = obj.tau_l{j};
                        case {'e', 'n'}
                            o = obj.tau_r{j};
                    end

                case 'tau_tot'
                    [e, tau] = obj.getBoundaryOperator({'e', 'tau'}, boundary);

                    I_dim = speye(obj.dim, obj.dim);
                    e_tot = kron(e, I_dim);
                    E = obj.E;
                    tau_tot = (e_tot'*E{1}*e*tau{1}')';
                    for i = 2:obj.dim
                        tau_tot = tau_tot + (e_tot'*E{i}*e*tau{i}')';
                    end
                    o = tau_tot;

                case 'H'
                    o = obj.H_boundary{j};

                case 'alpha'
                    % alpha = alpha(i,j) is the penalty strength for displacement BC.
                    e = obj.getBoundaryOperator('e', boundary);

                    LAMBDA = obj.LAMBDA;
                    MU = obj.MU;

                    dim = obj.dim;
                    theta_R = obj.theta_R{j};
                    theta_H = obj.theta_H{j};
                    theta_M = obj.theta_M{j};

                    a_lambda = dim/theta_H + 1/theta_R;
                    a_mu_i = 2/theta_M;
                    a_mu_ij = 2/theta_H + 1/theta_R;

                    d = @kroneckerDelta;  % Kronecker delta
                    db = @(i,j) 1-d(i,j); % Logical not of Kronecker delta
                    alpha = cell(obj.dim, obj.dim);

                    alpha_func = @(i,j) d(i,j)* a_lambda*LAMBDA ...
                                        + d(i,j)* a_mu_i*MU ...
                                        + db(i,j)*a_mu_ij*MU;
                    for i = 1:obj.dim
                        for l = 1:obj.dim
                            alpha{i,l} = d(i,l)*alpha_func(i,j)*e;
                        end
                    end

                    o = alpha;

                case 'alpha_tot'
                    % alpha = alpha(i,j) is the penalty strength for displacement BC.
                    [e, e_tot, alpha] = obj.getBoundaryOperator({'e', 'e_tot', 'alpha'}, boundary);
                    E = obj.E;
                    [m, n] = size(alpha{1,1});
                    alpha_tot = sparse(m*obj.dim, n*obj.dim);
                    for i = 1:obj.dim
                        for l = 1:obj.dim
                            alpha_tot = alpha_tot + (e_tot'*E{i}*e*alpha{i,l}'*E{l}')';
                        end
                    end
                    o = alpha_tot;
            end

        end

        % Returns square boundary quadrature matrix, of dimension
        % corresponding to the number of boundary points
        %
        % boundary -- string
        function H = getBoundaryQuadrature(obj, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            switch boundary
                case {'w','e'}
                    j = 1;
                case {'s','n'}
                    j = 2;
            end
            H = obj.H_boundary{j};
            I_dim = speye(obj.dim, obj.dim);
            H = kron(H, I_dim);
        end

        function N = size(obj)
            N = obj.dim*prod(obj.m);
        end
    end
end
