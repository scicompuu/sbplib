classdef Heat2dCurvilinear < scheme.Scheme

% Discretizes the Laplacian with variable coefficent, curvilinear,
% in the Heat equation way (i.e., the discretization matrix is not necessarily 
% symmetric)
% u_t = div * (kappa * grad u ) 
% opSet should be cell array of opSets, one per dimension. This
% is useful if we have periodic BC in one direction.

    properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing

        grid
        dim

        order % Order of accuracy for the approximation

        % Diagonal matrix for variable coefficients
        KAPPA % Variable coefficient

        D % Total operator
        D1 % First derivatives

        % Second derivatives
        D2_kappa

        H, Hi % Inner products
        e_l, e_r
        d1_l, d1_r % Normal derivatives at the boundary
        alpha % Vector of borrowing constants
        
        % Boundary inner products
        H_boundary_l, H_boundary_r 

        % Metric coefficients
        b % Cell matrix of size dim x dim
        J, Ji
        beta % Cell array of scale factors

        % Numerical boundary flux operators
        flux_l, flux_r

    end

    methods

        function obj = Heat2dCurvilinear(g ,order, kappa_fun, opSet)
            default_arg('opSet',{@sbp.D2Variable, @sbp.D2Variable});
            default_arg('kappa_fun', @(x,y) 0*x+1);
            dim = 2;

            kappa = grid.evalOn(g, kappa_fun);
            m = g.size();
            m_tot = g.N();

            % 1D operators
            ops = cell(dim,1);
            for i = 1:dim
                ops{i} = opSet{i}(m(i), {0, 1}, order);
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
            KAPPA = spdiag(kappa);
            obj.KAPPA = KAPPA;

            % Allocate
            obj.D1 = cell(dim,1);
            obj.D2_kappa = cell(dim,1);
            obj.e_l = cell(dim,1);
            obj.e_r = cell(dim,1);
            obj.d1_l = cell(dim,1);
            obj.d1_r = cell(dim,1);

            % D1
            obj.D1{1} = kron(D1{1},I{2});
            obj.D1{2} = kron(I{1},D1{2});

            % -- Metric coefficients ----
            coords = g.points();
            x = coords(:,1);
            y = coords(:,2);

            % Use non-periodic difference operators for metric even if opSet is periodic.
            xmax = max(ops{1}.x);
            ymax = max(ops{2}.x);
            opSetMetric{1} = sbp.D2Variable(m(1), {0, xmax}, order);
            opSetMetric{2} = sbp.D2Variable(m(2), {0, ymax}, order);
            D1Metric{1} = kron(opSetMetric{1}.D1, I{2});
            D1Metric{2} = kron(I{1}, opSetMetric{2}.D1); 

            x_xi = D1Metric{1}*x;
            x_eta = D1Metric{2}*x;
            y_xi = D1Metric{1}*y;
            y_eta = D1Metric{2}*y;

            J = x_xi.*y_eta - x_eta.*y_xi;

            b = cell(dim,dim);
            b{1,1} = y_eta./J;
            b{1,2} = -x_eta./J;
            b{2,1} = -y_xi./J;
            b{2,2} = x_xi./J;

            % Scale factors for boundary integrals
            beta = cell(dim,1);
            beta{1} = sqrt(x_eta.^2 + y_eta.^2);
            beta{2} = sqrt(x_xi.^2 + y_xi.^2);

            J = spdiag(J);
            Ji = inv(J);
            for i = 1:dim
                beta{i} = spdiag(beta{i});
                for j = 1:dim
                    b{i,j} = spdiag(b{i,j});
                end
            end
            obj.J = J;
            obj.Ji = Ji;
            obj.b = b;
            obj.beta = beta;
            %----------------------------

            % Boundary operators
            obj.e_l{1} = kron(e_l{1},I{2});
            obj.e_l{2} = kron(I{1},e_l{2});
            obj.e_r{1} = kron(e_r{1},I{2});
            obj.e_r{2} = kron(I{1},e_r{2});

            obj.d1_l{1} = kron(d1_l{1},I{2});
            obj.d1_l{2} = kron(I{1},d1_l{2});
            obj.d1_r{1} = kron(d1_r{1},I{2});
            obj.d1_r{2} = kron(I{1},d1_r{2});

            % D2 coefficients
            kappa_coeff = cell(dim,dim);
            for j = 1:dim
                obj.D2_kappa{j} = sparse(m_tot,m_tot); 
                kappa_coeff{j} = sparse(m_tot,1);
                for i = 1:dim
                    kappa_coeff{j} = kappa_coeff{j} + b{i,j}*J*b{i,j}*kappa;
                end
            end
            ind = grid.funcToMatrix(g, 1:m_tot);

            % x-dir
            j = 1;
            for col = 1:m(2)
                D_kappa = D2{1}(kappa_coeff{j}(ind(:,col)));

                p = ind(:,col);
                obj.D2_kappa{j}(p,p) = D_kappa;
            end

            % y-dir
            j = 2;
            for row = 1:m(1)
                D_kappa = D2{2}(kappa_coeff{j}(ind(row,:)));

                p = ind(row,:);
                obj.D2_kappa{j}(p,p) = D_kappa;
            end

            % Quadratures
            obj.H = kron(H{1},H{2});
            obj.Hi = inv(obj.H);
            obj.H_boundary_l = cell(dim,1);
            obj.H_boundary_l{1} = obj.e_l{1}'*beta{1}*obj.e_l{1}*H{2};
            obj.H_boundary_l{2} = obj.e_l{2}'*beta{2}*obj.e_l{2}*H{1};
            obj.H_boundary_r = cell(dim,1);
            obj.H_boundary_r{1} = obj.e_r{1}'*beta{1}*obj.e_r{1}*H{2};
            obj.H_boundary_r{2} = obj.e_r{2}'*beta{2}*obj.e_r{2}*H{1};

            %=== Differentiation matrix D (without SAT) ===
            D2_kappa = obj.D2_kappa;
            D1 = obj.D1;
            D = sparse(m_tot,m_tot);

            d = @kroneckerDelta;    % Kronecker delta
            db = @(i,j) 1-d(i,j); % Logical not of Kronecker delta

            % 2nd derivatives
            for j = 1:dim
                D = D + Ji*D2_kappa{j};
            end

            % Mixed terms
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        D = D + db(i,j)*Ji*D1{j}*b{i,j}*J*KAPPA*b{i,k}*D1{k};
                    end
                end
            end
            obj.D = D;
            %=========================================%

            % Normal flux operators for BC.
            flux_l = cell(dim,1);
            flux_r = cell(dim,1);

            d1_l = obj.d1_l;
            d1_r = obj.d1_r;
            e_l = obj.e_l;
            e_r = obj.e_r;

            % Loop over boundaries
            for j = 1:dim
                flux_l{j} = sparse(m_tot,m_tot);
                flux_r{j} = sparse(m_tot,m_tot);

                % Loop over dummy index
                for i = 1:dim
                    % Loop over dummy index
                    for k = 1:dim
                        flux_l{j} = flux_l{j} ...
                                  - beta{j}\b{i,j}*J*KAPPA*b{i,k}*( d(j,k)*e_l{k}*d1_l{k}' + db(j,k)*D1{k} );

                        flux_r{j} = flux_r{j} ...
                                  + beta{j}\b{i,j}*J*KAPPA*b{i,k}*( d(j,k)*e_r{k}*d1_r{k}' + db(j,k)*D1{k} );
                    end

                end
            end
            obj.flux_l = flux_l;
            obj.flux_r = flux_r;

            % Misc.
            obj.m = m;
            obj.h = g.scaling();
            obj.order = order;
            obj.grid = g;
            obj.dim = dim;
            obj.alpha = [ops{1}.borrowing.M.d1, ops{2}.borrowing.M.d1];

        end


        % Closure functions return the operators applied to the own domain to close the boundary
        % Penalty functions return the operators to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of boundary condition.
        %       data                is a function returning the data that should be applied at the boundary.
        %       neighbour_scheme    is an instance of Scheme that should be interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to interface to.
        function [closure, penalty] = boundary_condition(obj, boundary, type, symmetric, tuning)
            default_arg('type','Neumann');
            default_arg('symmetric', false);
            default_arg('tuning',1.2);

            % j is the coordinate direction of the boundary
            % nj: outward unit normal component. 
            % nj = -1 for west, south, bottom boundaries
            % nj = 1  for east, north, top boundaries
            [j, nj] = obj.get_boundary_number(boundary);
            switch nj
            case 1
                e = obj.e_r{j};
                flux = obj.flux_r{j};
                H_gamma = obj.H_boundary_r{j};
            case -1
                e = obj.e_l{j};
                flux = obj.flux_l{j};
                H_gamma = obj.H_boundary_l{j};
            end

            Hi = obj.Hi;
            Ji = obj.Ji;
            KAPPA = obj.KAPPA;
            kappa_gamma = e'*KAPPA*e; 
            h = obj.h(j);
            alpha = h*obj.alpha(j);

            switch type

            % Dirichlet boundary condition
            case {'D','d','dirichlet','Dirichlet'}

                if ~symmetric
                    closure = -Ji*Hi*flux'*e*H_gamma*(e' ); 
                    penalty = Ji*Hi*flux'*e*H_gamma;
                else
                    closure = Ji*Hi*flux'*e*H_gamma*(e' )...
                              -tuning*2/alpha*Ji*Hi*e*kappa_gamma*H_gamma*(e' ) ; 
                    penalty =  -Ji*Hi*flux'*e*H_gamma ...
                              +tuning*2/alpha*Ji*Hi*e*kappa_gamma*H_gamma;
                end

            % Normal flux boundary condition
            case {'N','n','neumann','Neumann'}
                    closure = -Ji*Hi*e*H_gamma*(e'*flux ); 
                    penalty =  Ji*Hi*e*H_gamma; 

            % Unknown boundary condition
            otherwise
                error('No such boundary condition: type = %s',type);
            end
        end

        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary)
            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            error('Interface not implemented');
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

        % Returns the coordinate number and outward normal component for the boundary specified by the string boundary.
        function [return_op] = get_boundary_operator(obj, op, boundary)

            switch boundary
                case {'w','W','west','West', 'e', 'E', 'east', 'East'}
                    j = 1;
                case {'s','S','south','South', 'n', 'N', 'north', 'North'}
                    j = 2;
                otherwise
                    error('No such boundary: boundary = %s',boundary);
            end

            switch op
                case 'e'
                    switch boundary
                        case {'w','W','west','West','s','S','south','South'}
                            return_op = obj.e_l{j};
                        case {'e', 'E', 'east', 'East','n', 'N', 'north', 'North'}
                            return_op = obj.e_r{j};
                    end
                case 'd'
                    switch boundary
                        case {'w','W','west','West','s','S','south','South'}
                            return_op = obj.d1_l{j};
                        case {'e', 'E', 'east', 'East','n', 'N', 'north', 'North'}
                            return_op = obj.d1_r{j};
                    end
                otherwise
                    error(['No such operator: operatr = ' op]);
            end

        end

        function N = size(obj)
            N = prod(obj.m);
        end
    end
end
