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

        H, Hi % Inner products
        phi % Borrowing constant for (d1 - e^T*D1) from R
        gamma % Borrowing constant for d1 from M
        H11 % First element of H
        e_l, e_r
        d1_l, d1_r % Normal derivatives at the boundary
        E % E{i}^T picks out component i
        
        H_boundary % Boundary inner products

        % Kroneckered norms and coefficients
        RHOi_kron
        Hi_kron
    end

    methods

        function obj = Elastic2dVariable(g ,order, lambda_fun, mu_fun, rho_fun, opSet)
            default_arg('opSet',{@sbp.D2Variable, @sbp.D2Variable});
            default_arg('lambda_fun', @(x,y) 0*x+1);
            default_arg('mu_fun', @(x,y) 0*x+1);
            default_arg('rho_fun', @(x,y) 0*x+1);
            dim = 2;

            assert(isa(g, 'grid.Cartesian'))

            lambda = grid.evalOn(g, lambda_fun);
            mu = grid.evalOn(g, mu_fun);
            rho = grid.evalOn(g, rho_fun);
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
                beta = ops{i}.borrowing.R.delta_D;
                obj.H11{i} = ops{i}.borrowing.H11;
                obj.phi{i} = beta/obj.H11{i};
                obj.gamma{i} = ops{i}.borrowing.M.d1;
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
            %=========================================%

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

                % Loop over components
                for i = 1:dim
                    tau_l{j}{i} = sparse(m_tot,dim*m_tot);
                    tau_r{j}{i} = sparse(m_tot,dim*m_tot);
                    for k = 1:dim
                        T_l{j}{i,k} = ... 
                        -d(i,j)*LAMBDA*(d(i,k)*e_l{k}*d1_l{k}' + db(i,k)*D1{k})...
                        -d(j,k)*MU*(d(i,j)*e_l{i}*d1_l{i}' + db(i,j)*D1{i})... 
                        -d(i,k)*MU*e_l{j}*d1_l{j}';

                        T_r{j}{i,k} = ... 
                        d(i,j)*LAMBDA*(d(i,k)*e_r{k}*d1_r{k}' + db(i,k)*D1{k})...
                        +d(j,k)*MU*(d(i,j)*e_r{i}*d1_r{i}' + db(i,j)*D1{i})... 
                        +d(i,k)*MU*e_r{j}*d1_r{j}';

                        tau_l{j}{i} = tau_l{j}{i} + T_l{j}{i,k}*E{k}';
                        tau_r{j}{i} = tau_r{j}{i} + T_r{j}{i,k}*E{k}';
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

        end


        % Closure functions return the operators applied to the own domain to close the boundary
        % Penalty functions return the operators to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a cell array of strings specifying the type of boundary condition for each component.
        %       data                is a function returning the data that should be applied at the boundary.
        %       neighbour_scheme    is an instance of Scheme that should be interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to interface to.
        function [closure, penalty] = boundary_condition(obj, boundary, type, parameter)
            default_arg('type',{'free','free'});
            default_arg('parameter', []);

            if ~iscell(type)
                type = {type, type};
            end

            % j is the coordinate direction of the boundary
            j = obj.get_boundary_number(boundary);
            [e, T, tau, H_gamma] = obj.get_boundary_operator({'e','T','tau','H'}, boundary);

            E = obj.E;
            Hi = obj.Hi;
            LAMBDA = obj.LAMBDA;
            MU = obj.MU;
            RHOi = obj.RHOi;

            dim = obj.dim;
            m_tot = obj.grid.N();

            % Preallocate
            closure = sparse(dim*m_tot, dim*m_tot);
            penalty = cell(dim,1);
            for k = 1:dim
                penalty{k} = sparse(dim*m_tot, m_tot/obj.m(j));
            end

            % Loop over components that we (potentially) have different BC on
            for k = 1:dim
                switch type{k}

                % Dirichlet boundary condition
                case {'D','d','dirichlet','Dirichlet'}

                    tuning = 1.2;
                    phi = obj.phi{j};
                    h = obj.h(j);
                    h11 = obj.H11{j}*h;
                    gamma = obj.gamma{j};

                    a_lambda = dim/h11 + 1/(h11*phi);
                    a_mu_i = 2/(gamma*h);
                    a_mu_ij = 2/h11 + 1/(h11*phi);

                    d = @kroneckerDelta;  % Kronecker delta
                    db = @(i,j) 1-d(i,j); % Logical not of Kronecker delta
                    alpha = @(i,j) tuning*( d(i,j)* a_lambda*LAMBDA ...
                                          + d(i,j)* a_mu_i*MU ...
                                          + db(i,j)*a_mu_ij*MU ); 

                    % Loop over components that Dirichlet penalties end up on
                    for i = 1:dim
                        C = T{k,i};
                        A = -d(i,k)*alpha(i,j);
                        B = A + C;
                        closure = closure + E{i}*RHOi*Hi*B'*e*H_gamma*(e'*E{k}' ); 
                        penalty{k} = penalty{k} - E{i}*RHOi*Hi*B'*e*H_gamma;
                    end 

                % Free boundary condition
                case {'F','f','Free','free','traction','Traction','t','T'}
                        closure = closure - E{k}*RHOi*Hi*e*H_gamma* (e'*tau{k} ); 
                        penalty{k} = penalty{k} + E{k}*RHOi*Hi*e*H_gamma;

                % Unknown boundary condition
                otherwise
                    error('No such boundary condition: type = %s',type);
                end
            end
        end

        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary)
            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            % Operators without subscripts are from the own domain.
            tuning = 1.2;

            % j is the coordinate direction of the boundary
            j = obj.get_boundary_number(boundary);
            j_v = neighbour_scheme.get_boundary_number(neighbour_boundary);

            % Get boundary operators
            [e, T, tau, H_gamma] = obj.get_boundary_operator({'e','T','tau','H'}, boundary);
            [e_v, tau_v] = neighbour_scheme.get_boundary_operator({'e','tau'}, neighbour_boundary);

            % Operators and quantities that correspond to the own domain only
            Hi = obj.Hi;
            RHOi = obj.RHOi;
            dim = obj.dim;
        
            %--- Other operators ----
            m_tot_u = obj.grid.N();
            E = obj.E;
            LAMBDA_u = obj.LAMBDA;
            MU_u = obj.MU;
            lambda_u = e'*LAMBDA_u*e;
            mu_u = e'*MU_u*e;

            m_tot_v = neighbour_scheme.grid.N();
            E_v = neighbour_scheme.E;
            LAMBDA_v = neighbour_scheme.LAMBDA;
            MU_v = neighbour_scheme.MU;
            lambda_v = e_v'*LAMBDA_v*e_v;
            mu_v = e_v'*MU_v*e_v;
            %-------------------------
            
            % Borrowing constants
            phi_u = obj.phi{j};
            h_u = obj.h(j);
            h11_u = obj.H11{j}*h_u;
            gamma_u = obj.gamma{j};

            phi_v = neighbour_scheme.phi{j_v};
            h_v = neighbour_scheme.h(j_v);
            h11_v = neighbour_scheme.H11{j_v}*h_v;
            gamma_v = neighbour_scheme.gamma{j_v};

            % E > sum_i 1/(2*alpha_ij)*(tau_i)^2
            function [alpha_ii, alpha_ij] = computeAlpha(phi,h,h11,gamma,lambda,mu) 
                th1 = h11/(2*dim);
                th2 = h11*phi/2;
                th3 = h*gamma;
                a1 = ( (th1 + th2)*th3*lambda + 4*th1*th2*mu ) / (2*th1*th2*th3);
                a2 = ( 16*(th1 + th2)*lambda*mu ) / (th1*th2*th3);
                alpha_ii = a1 + sqrt(a2 + a1^2);

                alpha_ij = 2/h11 + 1/(phi*h11);
            end

            [alpha_ii_u, alpha_ij_u] = computeAlpha(phi_u,h_u,h11_u,gamma_u,lambda_u,mu_u);
            [alpha_ii_v, alpha_ij_v] = computeAlpha(phi_v,h_v,h11_v,gamma_v,lambda_v,mu_v);  
            sigma_ii = tuning*(alpha_ii_u + alpha_ii_v)/4;
            sigma_ij = tuning*(alpha_ij_u + alpha_ij_v)/4;

            d = @kroneckerDelta;  % Kronecker delta
            db = @(i,j) 1-d(i,j); % Logical not of Kronecker delta
            sigma = @(i,j) tuning*(d(i,j)*sigma_ii + db(i,j)*sigma_ij);

            % Preallocate
            closure = sparse(dim*m_tot_u, dim*m_tot_u);
            penalty = sparse(dim*m_tot_u, dim*m_tot_v);

            % Loop over components that penalties end up on
            for i = 1:dim
                closure = closure - E{i}*RHOi*Hi*e*sigma(i,j)*H_gamma*e'*E{i}';
                penalty = penalty + E{i}*RHOi*Hi*e*sigma(i,j)*H_gamma*e_v'*E_v{i}';

                closure = closure - 1/2*E{i}*RHOi*Hi*e*H_gamma*e'*tau{i};
                penalty = penalty - 1/2*E{i}*RHOi*Hi*e*H_gamma*e_v'*tau_v{i};

                % Loop over components that we have interface conditions on
                for k = 1:dim
                    closure = closure + 1/2*E{i}*RHOi*Hi*T{k,i}'*e*H_gamma*e'*E{k}'; 
                    penalty = penalty - 1/2*E{i}*RHOi*Hi*T{k,i}'*e*H_gamma*e_v'*E_v{k}'; 
                end 
            end
        end

        % Returns the coordinate number and outward normal component for the boundary specified by the string boundary.
        function [j, nj] = get_boundary_number(obj, boundary)

            switch boundary
                case {'w','W','west','West', 'e', 'E', 'east', 'East'}
                    j = 1;
                case {'s','S','south','South', 'n', 'N', 'north', 'North'}
                    j = 2;
                otherwise
                    error('No such boundary: boundary = %s',boundary);
            end

            switch boundary
                case {'w','W','west','West','s','S','south','South'}
                    nj = -1;
                case {'e', 'E', 'east', 'East','n', 'N', 'north', 'North'}
                    nj = 1;
            end
        end

        % Returns the boundary operator op for the boundary specified by the string boundary.
        % op: may be a cell array of strings
        function [varargout] = get_boundary_operator(obj, op, boundary)

            switch boundary
                case {'w','W','west','West', 'e', 'E', 'east', 'East'}
                    j = 1;
                case {'s','S','south','South', 'n', 'N', 'north', 'North'}
                    j = 2;
                otherwise
                    error('No such boundary: boundary = %s',boundary);
            end

            if ~iscell(op)
                op = {op};
            end

            for i = 1:length(op)
                switch op{i}
                    case 'e'
                        switch boundary
                            case {'w','W','west','West','s','S','south','South'}
                                varargout{i} = obj.e_l{j};
                            case {'e', 'E', 'east', 'East','n', 'N', 'north', 'North'}
                                varargout{i} = obj.e_r{j};
                        end
                    case 'd'
                        switch boundary
                            case {'w','W','west','West','s','S','south','South'}
                                varargout{i} = obj.d1_l{j};
                            case {'e', 'E', 'east', 'East','n', 'N', 'north', 'North'}
                                varargout{i} = obj.d1_r{j};
                        end
                    case 'H'
                        varargout{i} = obj.H_boundary{j};
                    case 'T'
                        switch boundary
                            case {'w','W','west','West','s','S','south','South'}
                                varargout{i} = obj.T_l{j};
                            case {'e', 'E', 'east', 'East','n', 'N', 'north', 'North'}
                                varargout{i} = obj.T_r{j};
                        end
                    case 'tau'
                        switch boundary
                            case {'w','W','west','West','s','S','south','South'}
                                varargout{i} = obj.tau_l{j};
                            case {'e', 'E', 'east', 'East','n', 'N', 'north', 'North'}
                                varargout{i} = obj.tau_r{j};
                        end                        
                    otherwise
                        error(['No such operator: operator = ' op{i}]);
                end
            end

        end

        function N = size(obj)
            N = prod(obj.m);
        end
    end
end
