classdef elasticShearVariable < scheme.Scheme
    properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing

        grid

        order % Order accuracy for the approximation

        a % Variable coefficient lambda of the operator
        rho % Density

        D % Total operator
        D1 % First derivatives
        D2 % Second derivatives
        H, Hi % Inner products
        e_l, e_r
        d_l, d_r % Normal derivatives at the boundary
        H_boundary % Boundary inner products
    end

    methods
        % Implements the shear part of the elastic wave equation, i.e.
        % rho u_{i,tt} = d_i a d_j u_j + d_j a d_j u_i
        % where a = lambda.

        function obj = elasticShearVariable(g ,order, a_fun, rho_fun, opSet)
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

            % 1D operators
            ops = cell(dim,1);
            for i = 1:dim
                ops{i} = opSet(m(i), {0, 1}, order);
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
                I{i} = speye{m(i));
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
            RHO = spdiag(rho);

            obj.D1 = cell(dim,1);
            obj.D2 = cell(dim,1);
            obj.e_l = cell(dim,1);
            obj.e_r = cell(dim,1);
            obj.d1_l = cell(dim,1);
            obj.d1_r = cell(dim,1);

            % D1
            obj.D1{1} = kron{D1{1},I(2)};
            obj.D1{2} = kron{I{1},D1(2)};

            % Boundary operators
            obj.e_l{1} = kron{e_l{1},I(2)};
            obj.e_l{2} = kron{I{1},e_l(2)};
            obj.e_r{1} = kron{e_r{1},I(2)};
            obj.e_r{2} = kron{I{1},e_r(2)};

            obj.d1_l{1} = kron{d1_l{1},I(2)};
            obj.d1_l{2} = kron{I{1},d1_l(2)};
            obj.d1_r{1} = kron{d1_r{1},I(2)};
            obj.d1_r{2} = kron{I{1},d1_r(2)};

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
            obj.H_boundary = cell(dim,1);
            obj.H_boundary{1} = H{2};
            obj.H_boundary{2} = H{1};

            % Boundary coefficient matrices and quadratures
            obj.A_boundary_l = cell(dim,1);
            obj.A_boundary_r = cell(dim,1);
            for i = 1:dim
                obj.A_boundary_l{i} = e_l{i}'*A*e_l{i};
                obj.A_boundary_r{i} = e_r{i}'*A*e_r{i};
            end

            % E{i}^T picks out component i.
            E = cell(dim,1);
            I = speye{mtot,mtot};
            for i = 1:dim
                e = sparse(dim,1);
                e(i) = 1;
                E{i} = kron(e,I)
            end
            obj.E = E;

            % Differentiation matrix D (without SAT)
            D = 0;
            d = @kroneckerDelta;    % Kronecker delta
            db = @(i,j) 1-dij(i,j); % Logical not of Kronecker delta
            for i = 1:dim
                for j = 1:dim
                    D = D + E{i}*inv(rho)*( d(i,j)*D2{i}*E{j}' +...
                                            db(i,j)*D1{j}*A*D1{i}*E{j}' + ...
                                            D2{j}*E{i}' ...
                                          );
                end
            end
            obj.D = D;
            %=========================================%

            

            % Misc.
            obj.m = m;
            obj.h = h;
            obj.order = order;
            obj.grid = g;

            obj.a = a;
            obj.b = b;

            % obj.gamm_u = h_u*ops_u.borrowing.M.d1;
            % obj.gamm_v = h_v*ops_v.borrowing.M.d1;
        end


        % Closure functions return the operators applied to the own domain to close the boundary
        % Penalty functions return the operators to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of boundary condition if there are several.
        %       data                is a function returning the data that should be applied at the boundary.
        %       neighbour_scheme    is an instance of Scheme that should be interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to interface to.
        function [closure, penalty] = boundary_condition(obj, boundary, type, parameter)
            default_arg('type','free');
            default_arg('parameter', []);

            delta = @kroneckerDelta;    % Kronecker delta
            delta_b = @(i,j) 1-dij(i,j); % Logical not of Kronecker delta

            % j is the coordinate direction of the boundary
            % nj: outward unit normal component. 
            % nj = -1 for west, south, bottom boundaries
            % nj = 1  for east, north, top boundaries
            [j, nj] = obj.get_boundary_number(boundary);
            switch nj
            case 1
                e = obj.e_r;
                d = obj.d_r;
            case -1
                e = obj.e_l;
                d = obj.d_l;
            end

            E = obj.E;
            Hi = obj.Hi;
            H_gamma = obj.H_boundary{j};
            A = obj.A;

            switch type
                % Dirichlet boundary condition
                case {'D','d','dirichlet'}
                    error('Dirichlet not implemented')
                    tuning = 1.2;
                    % tuning = 20.2;

                    b1 = gamm*obj.lambda./obj.a11.^2;
                    b2 = gamm*obj.lambda./obj.a22.^2;

                    tau1 = tuning * spdiag(-1./b1 - 1./b2);
                    tau2 = 1;

                    tau = (tau1*e + tau2*d)*H_b;

                    closure =  obj.a*obj.Hi*tau*e';
                    penalty = -obj.a*obj.Hi*tau;


                % Free boundary condition
                case {'F','f','Free','free'}
                    closure = 0;
                    penalty = 0;
                    % Loop over components
                    for i = 1:3
                        closure = closure + E{i}*(-sign)*Hi*e{j}*H_gamma*(...
                                e{j}'*A*e{j}*d{j}'*E{i}' + ...
                                delta(i,j)*e{j}'*A*e{i}*d{i}*E{j}' + ...
                                delta_b(i,j)*A*D1{i}*E{j}' ...
                                );
                        penalty = penalty - E{i}*(-sign)*Hi*e{j}*H_gamma;
                    end

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

        % Ruturns the coordinate number and outward normal component for the boundary specified by the string boundary.
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

        function N = size(obj)
            N = prod(obj.m);
        end
    end
end
