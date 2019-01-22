classdef Heat2dVariable < scheme.Scheme

% Discretizes the Laplacian with variable coefficent,
% In the Heat equation way (i.e., the discretization matrix is not necessarily
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

        H_boundary % Boundary inner products

    end

    methods

        function obj = Heat2dVariable(g ,order, kappa_fun, opSet)
            default_arg('opSet',{@sbp.D2Variable, @sbp.D2Variable});
            default_arg('kappa_fun', @(x,y) 0*x+1);
            dim = 2;

            assert(isa(g, 'grid.Cartesian'))

            kappa = grid.evalOn(g, kappa_fun);
            m = g.size();
            m_tot = g.N();

            h = g.scaling();
            lim = g.lim;

            % 1D operators
            ops = cell(dim,1);
            for i = 1:dim
                ops{i} = opSet{i}(m(i), lim{i}, order);
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

            obj.D1 = cell(dim,1);
            obj.D2_kappa = cell(dim,1);
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
                obj.D2_kappa{i} = sparse(m_tot);
            end
            ind = grid.funcToMatrix(g, 1:m_tot);

            for i = 1:m(2)
                D_kappa = D2{1}(kappa(ind(:,i)));
                p = ind(:,i);
                obj.D2_kappa{1}(p,p) = D_kappa;
            end

            for i = 1:m(1)
                D_kappa = D2{2}(kappa(ind(i,:)));
                p = ind(i,:);
                obj.D2_kappa{2}(p,p) = D_kappa;
            end

            % Quadratures
            obj.H = kron(H{1},H{2});
            obj.Hi = inv(obj.H);
            obj.H_boundary = cell(dim,1);
            obj.H_boundary{1} = H{2};
            obj.H_boundary{2} = H{1};

            % Differentiation matrix D (without SAT)
            D2_kappa = obj.D2_kappa;
            D1 = obj.D1;
            D = sparse(m_tot,m_tot);
            for i = 1:dim
                D = D + D2_kappa{i};
            end
            obj.D = D;
            %=========================================%

            % Misc.
            obj.m = m;
            obj.h = h;
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

            % nj: outward unit normal component.
            % nj = -1 for west, south, bottom boundaries
            % nj = 1  for east, north, top boundaries
            nj = obj.getBoundarySign(boundary);

            Hi = obj.Hi;
            [e, d] = obj.getBoundaryOperator({'e', 'd'}, boundary);
            H_gamma = obj.getBoundaryQuadrature(boundary);
            alpha = obj.getBoundaryBorrowing(boundary);

            KAPPA = obj.KAPPA;
            kappa_gamma = e'*KAPPA*e;

            switch type

            % Dirichlet boundary condition
            case {'D','d','dirichlet','Dirichlet'}

                if ~symmetric
                    closure = -nj*Hi*d*kappa_gamma*H_gamma*(e' );
                    penalty =  nj*Hi*d*kappa_gamma*H_gamma;
                else
                    closure = nj*Hi*d*kappa_gamma*H_gamma*(e' )...
                              -tuning*2/alpha*Hi*e*kappa_gamma*H_gamma*(e' ) ;
                    penalty =  -nj*Hi*d*kappa_gamma*H_gamma ...
                              +tuning*2/alpha*Hi*e*kappa_gamma*H_gamma;
                end

            % Free boundary condition
            case {'N','n','neumann','Neumann'}
                    closure = -nj*Hi*e*kappa_gamma*H_gamma*(d' );
                    penalty =  Hi*e*kappa_gamma*H_gamma;
                    % penalty is for normal derivative and not for derivative, hence the sign.

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

        % Returns the boundary operator op for the boundary specified by the string boundary.
        % op        -- string or a cell array of strings
        % boundary  -- string
        function varargout = getBoundaryOperator(obj, op, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            if ~iscell(op)
                op = {op};
            end

            for i = 1:numel(op)
                switch op{i}
                case 'e'
                    switch boundary
                    case 'w'
                        e = obj.e_l{1};
                    case 'e'
                        e = obj.e_r{1};
                    case 's'
                        e = obj.e_l{2};
                    case 'n'
                        e = obj.e_r{2};
                    end
                    varargout{i} = e;

                case 'd'
                    switch boundary
                    case 'w'
                        d = obj.d1_l{1};
                    case 'e'
                        d = obj.d1_r{1};
                    case 's'
                        d = obj.d1_l{2};
                    case 'n'
                        d = obj.d1_r{2};
                    end
                    varargout{i} = d;
                end
            end
        end

        % Returns square boundary quadrature matrix, of dimension
        % corresponding to the number of boundary points
        %
        % boundary -- string
        function H_b = getBoundaryQuadrature(obj, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            switch boundary
                case 'w'
                    H_b = obj.H_boundary{1};
                case 'e'
                    H_b = obj.H_boundary{1};
                case 's'
                    H_b = obj.H_boundary{2};
                case 'n'
                    H_b = obj.H_boundary{2};
            end
        end

        % Returns the boundary sign. The right boundary is considered the positive boundary
        % boundary -- string
        function s = getBoundarySign(obj, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            switch boundary
                case {'e','n'}
                    s = 1;
                case {'w','s'}
                    s = -1;
            end
        end

        % Returns borrowing constant gamma*h
        % boundary -- string
        function gamm = getBoundaryBorrowing(obj, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            switch boundary
                case {'w','e'}
                    gamm = obj.h(1)*obj.alpha(1);
                case {'s','n'}
                    gamm = obj.h(2)*obj.alpha(2);
            end
        end

        function N = size(obj)
            N = prod(obj.m);
        end
    end
end
