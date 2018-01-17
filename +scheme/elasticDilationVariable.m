classdef elasticDilationVariable < scheme.Scheme
    properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing

        grid
        dim

        order % Order accuracy for the approximation

        A % Variable coefficient lambda of the operator (as diagonal matrix here)
        RHO % Density (as diagonal matrix here)

        D % Total operator
        D1 % First derivatives
        D2 % Second derivatives
        Div % Divergence operator used for BC
        H, Hi % Inner products
        phi % Borrowing constant for (d1 - e^T*D1) from R
        H11 % First element of H
        e_l, e_r
        d1_l, d1_r % Normal derivatives at the boundary
        E % E{i}^T picks out component i
        
        H_boundary % Boundary inner products

        A_boundary_l % Variable coefficient at boundaries
        A_boundary_r % 
    end

    methods
        % Implements the shear part of the elastic wave equation, i.e.
        % rho u_{i,tt} = d_i a d_j u_j + d_j a d_j u_i
        % where a = mu.

        function obj = elasticDilationVariable(g ,order, a_fun, rho_fun, opSet)
            default_arg('opSet',@sbp.D2Variable);
            default_arg('a_fun', @(x,y) 0*x+1);
            default_arg('rho_fun', @(x,y) 0*x+1);
            dim = 2;

            assert(isa(g, 'grid.Cartesian'))

            a = grid.evalOn(g, a_fun);
            rho = grid.evalOn(g, rho_fun);
            m = g.size();
            m_tot = g.N();

            h = g.scaling();
            L = (m-1).*h;

            % 1D operators
            ops = cell(dim,1);
            for i = 1:dim
                ops{i} = opSet(m(i), {0, L(i)}, order);
            end

            % Borrowing constants
            beta = ops{1}.borrowing.R.delta_D;
            obj.H11 = ops{1}.borrowing.H11;
            obj.phi = beta/obj.H11;

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
            A = spdiag(a);
            obj.A = A;
            RHO = spdiag(rho);
            obj.RHO = RHO;


            obj.D1 = cell(dim,1);
            obj.D2 = cell(dim,1);
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
                obj.D2{i} = sparse(m_tot);
            end
            ind = grid.funcToMatrix(g, 1:m_tot);

            for i = 1:m(2)
                D = D2{1}(a(ind(:,i)));
                p = ind(:,i);
                obj.D2{1}(p,p) = D;
            end

            for i = 1:m(1)
                D = D2{2}(a(ind(i,:)));
                p = ind(i,:);
                obj.D2{2}(p,p) = D;
            end

            % Quadratures
            obj.H = kron(H{1},H{2});
            obj.Hi = inv(obj.H);
            obj.H_boundary = cell(dim,1);
            obj.H_boundary{1} = H{2};
            obj.H_boundary{2} = H{1};

            % Boundary coefficient matrices and quadratures
            obj.A_boundary_l = cell(dim,1);
            obj.A_boundary_r = cell(dim,1);
            for i = 1:dim
                obj.A_boundary_l{i} = obj.e_l{i}'*A*obj.e_l{i};
                obj.A_boundary_r{i} = obj.e_r{i}'*A*obj.e_r{i};
            end

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
            D2 = obj.D2;
            D1 = obj.D1;
            D = sparse(dim*m_tot,dim*m_tot);
            d = @kroneckerDelta;    % Kronecker delta
            db = @(i,j) 1-d(i,j); % Logical not of Kronecker delta
            for i = 1:dim
                for j = 1:dim
                    D = D + E{i}*inv(RHO)*( d(i,j)*D2{i}*E{j}' +...
                                            db(i,j)*D1{i}*A*D1{j}*E{j}' ...
                                          );
                end
            end
            obj.D = D;
            %=========================================%

            % Divergence operator for BC
            Div = cell(dim,1);
            for i = 1:dim
                Div{i} = sparse(m_tot,dim*m_tot);
                for j = 1:dim
                    Div{i} = Div{i} + d(i,j)*(obj.e_l{i}*obj.d1_l{i}' + obj.e_r{i}*obj.d1_r{i}')*E{j}' ...
                              + db(i,j)*obj.D1{j}*E{j}';
                end
            end
            obj.Div = Div;

            % Misc.
            obj.m = m;
            obj.h = h;
            obj.order = order;
            obj.grid = g;
            obj.dim = dim;

        end


        % Closure functions return the operators applied to the own domain to close the boundary
        % Penalty functions return the operators to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        % Here penalty{i,j} enforces data component j on solution component i
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of boundary condition if there are several.
        %       data                is a function returning the data that should be applied at the boundary.
        %       neighbour_scheme    is an instance of Scheme that should be interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to interface to.
        function [closure, penalty] = boundary_condition(obj, boundary, type, parameter)
            default_arg('type','free');
            default_arg('parameter', []);

            delta = @kroneckerDelta;    % Kronecker delta
            delta_b = @(i,j) 1-delta(i,j); % Logical not of Kronecker delta

            % j is the coordinate direction of the boundary
            % nj: outward unit normal component. 
            % nj = -1 for west, south, bottom boundaries
            % nj = 1  for east, north, top boundaries
            [j, nj] = obj.get_boundary_number(boundary);
            switch nj
            case 1
                e = obj.e_r;
                d = obj.d1_r;
            case -1
                e = obj.e_l;
                d = obj.d1_l;
            end

            E = obj.E;
            Hi = obj.Hi;
            H_gamma = obj.H_boundary{j};
            A = obj.A;
            RHO = obj.RHO;
            D1 = obj.D1;

            phi = obj.phi;
            H11 = obj.H11;
            h = obj.h;

            switch type
                % Dirichlet boundary condition
                case {'D','d','dirichlet'}
                    error('Not implemented');
                    tuning = 1.2;
                    phi = obj.phi;

                    closures = cell(obj.dim,1);
                    penalties = cell(obj.dim,obj.dim);
                    % Loop over components
                    for i = 1:obj.dim
                        H_gamma_i = obj.H_boundary{i};
                        sigma_ij = tuning*delta(i,j)*2/(gamma*h(j)) +...
                                   tuning*delta_b(i,j)*(2/(H11*h(j)) + 2/(H11*h(j)*phi));

                        ci = E{i}*inv(RHO)*nj*Hi*...
                                ( (e{j}*H_gamma*e{j}'*A*e{j}*d{j}')'*E{i}' + ...
                                  delta(i,j)*(e{j}*H_gamma*e{j}'*A*e{j}*d{j}')'*E{j}' ...
                                   ) ...
                                - sigma_ij*E{i}*inv(RHO)*Hi*A*e{j}*H_gamma*e{j}'*E{i}';

                        cj = E{j}*inv(RHO)*nj*Hi*...
                                ( delta_b(i,j)*(e{j}*H_gamma*e{j}'*A*D1{i})'*E{i}' ...
                                   );

                        if isempty(closures{i})
                            closures{i} = ci;
                        else
                            closures{i} = closures{i} + ci;
                        end

                        if isempty(closures{j})
                            closures{j} = cj;
                        else
                            closures{j} = closures{j} + cj;
                        end
   
                        pii = - E{i}*inv(RHO)*nj*Hi*...
                                ( (H_gamma*e{j}'*A*e{j}*d{j}')' + ...
                                  delta(i,j)*(H_gamma*e{j}'*A*e{j}*d{j}')' ...
                                   ) ...
                                + sigma_ij*E{i}*inv(RHO)*Hi*A*e{j}*H_gamma;

                        pji = - E{j}*inv(RHO)*nj*Hi*...
                                ( delta_b(i,j)*(H_gamma*e{j}'*A*D1{i})' );

                        % Dummies
                        pij = - 0*E{i}*e{j};
                        pjj = - 0*E{j}*e{j};

                        if isempty(penalties{i,i})
                            penalties{i,i} = pii;
                        else
                            penalties{i,i} = penalties{i,i} + pii;
                        end

                        if isempty(penalties{j,i})
                            penalties{j,i} = pji;
                        else
                            penalties{j,i} = penalties{j,i} + pji;
                        end

                        if isempty(penalties{i,j})
                            penalties{i,j} = pij;
                        else
                            penalties{i,j} = penalties{i,j} + pij;
                        end

                        if isempty(penalties{j,j})
                            penalties{j,j} = pjj;
                        else
                            penalties{j,j} = penalties{j,j} + pjj;
                        end
                    end
                    [rows, cols] = size(closures{1});
                    closure = sparse(rows, cols);
                    for i = 1:obj.dim
                        closure = closure + closures{i};
                    end
                    penalty = penalties;

                % Free boundary condition
                case {'F','f','Free','free'}
                    closures = cell(obj.dim,1);
                    penalties = cell(obj.dim,obj.dim);

                    % Divergence operator
                    Div = obj.Div{j};

                    % Loop over components
                    %for i = 1:obj.dim
                        closure = -nj*E{j}*inv(RHO)*Hi*e{j} ...
                                     *H_gamma*e{j}'*A*e{j}*e{j}'*Div;
                        penalty = nj*E{j}*inv(RHO)*Hi*e{j} ...
                                     *H_gamma*e{j}'*A*e{j};
                    %end
                    % [rows, cols] = size(closures{1});
                    % closure = sparse(rows, cols);
                    % for i = 1:obj.dim
                    %     closure = closure + closures{i};
                    %     for j = 1:obj.dim
                    %         if i~=j
                    %             [rows cols] = size(penalties{j,j});
                    %             penalties{i,j} = sparse(rows,cols);
                    %         end
                    %     end
                    % end
                    % penalty = penalties;


                % Unknown boundary condition
                otherwise
                    error('No such boundary condition: type = %s',type);
            end
        end

        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary)
            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            tuning = 1.2;
            % tuning = 20.2;
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
                            return_op = obj.d_l{j};
                        case {'e', 'E', 'east', 'East','n', 'N', 'north', 'North'}
                            return_op = obj.d_r{j};
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
