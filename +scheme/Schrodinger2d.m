classdef Schrodinger2d < scheme.Scheme

% Discretizes the Laplacian with constant coefficent,
% in the Schrödinger equation way (i.e., the discretization matrix is not necessarily
% definite)
% u_t = a*i*Laplace u
% opSet should be cell array of opSets, one per dimension. This
% is useful if we have periodic BC in one direction.

    properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing

        grid
        dim

        order % Order of accuracy for the approximation

        % Diagonal matrix for variable coefficients
        a % Constant coefficient

        D % Total operator
        D1 % First derivatives

        % Second derivatives
        D2

        H, Hi % Inner products
        e_l, e_r
        d1_l, d1_r % Normal derivatives at the boundary
        e_w, e_e, e_s, e_n
        d_w, d_e, d_s, d_n

        H_boundary % Boundary inner products

    end

    methods

        function obj = Schrodinger2d(g ,order, a, opSet)
            default_arg('opSet',{@sbp.D2Variable, @sbp.D2Variable});
            default_arg('a',1);
            dim = 2;

            assert(isa(g, 'grid.Cartesian'))
            if isa(a, 'function_handle')
                a = grid.evalOn(g, a);
                a = spdiag(a);
            end

            m = g.size();
            m_tot = g.N();

            h = g.scaling();
            xlim = {g.x{1}(1), g.x{1}(end)};
            ylim = {g.x{2}(1), g.x{2}(end)};
            lim = {xlim, ylim};

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

            % Constant coeff D2
            for i = 1:dim
                D2{i} = D2{i}(ones(m(i),1));
            end

            %====== Assemble full operators ========
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
            obj.D2{1} = kron(D2{1},I{2});
            obj.D2{2} = kron(I{1},D2{2});

            % Quadratures
            obj.H = kron(H{1},H{2});
            obj.Hi = inv(obj.H);
            obj.H_boundary = cell(dim,1);
            obj.H_boundary{1} = H{2};
            obj.H_boundary{2} = H{1};

            % Differentiation matrix D (without SAT)
            D2 = obj.D2;
            D = sparse(m_tot,m_tot);
            for j = 1:dim
                D = D + a*1i*D2{j};
            end
            obj.D = D;
            %=========================================%

            % Misc.
            obj.m = m;
            obj.h = h;
            obj.order = order;
            obj.grid = g;
            obj.dim = dim;
            obj.a = a;
            obj.e_w = obj.e_l{1};
            obj.e_e = obj.e_r{1};
            obj.e_s = obj.e_l{2};
            obj.e_n = obj.e_r{2};
            obj.d_w = obj.d1_l{1};
            obj.d_e = obj.d1_r{1};
            obj.d_s = obj.d1_l{2};
            obj.d_n = obj.d1_r{2};

        end


        % Closure functions return the operators applied to the own domain to close the boundary
        % Penalty functions return the operators to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of boundary condition.
        %       data                is a function returning the data that should be applied at the boundary.
        %       neighbour_scheme    is an instance of Scheme that should be interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to interface to.
        function [closure, penalty] = boundary_condition(obj, boundary, type, parameter)
            default_arg('type','Neumann');
            default_arg('parameter', []);

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

            Hi = obj.Hi;
            H_gamma = obj.H_boundary{j};
            a = e{j}'*obj.a*e{j};

            switch type

            % Dirichlet boundary condition
            case {'D','d','dirichlet','Dirichlet'}
                    closure =  nj*Hi*d{j}*a*1i*H_gamma*(e{j}' );
                    penalty = -nj*Hi*d{j}*a*1i*H_gamma;

            % Free boundary condition
            case {'N','n','neumann','Neumann'}
                    closure = -nj*Hi*e{j}*a*1i*H_gamma*(d{j}' );
                    penalty =  nj*Hi*e{j}*a*1i*H_gamma;

            % Unknown boundary condition
            otherwise
                error('No such boundary condition: type = %s',type);
            end
        end

        % type     Struct that specifies the interface coupling.
        %          Fields:
        %          -- interpolation:    type of interpolation, default 'none'
        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary,type)

            defaultType.interpolation = 'none';
            default_struct('type', defaultType);

            switch type.interpolation
            case {'none', ''}
                [closure, penalty] = interfaceStandard(obj,boundary,neighbour_scheme,neighbour_boundary);
            case {'op','OP'}
                [closure, penalty] = interfaceNonConforming(obj,boundary,neighbour_scheme,neighbour_boundary,type);
            otherwise
                error('Unknown type of interpolation: %s ', type.interpolation);
            end
        end

        function [closure, penalty] = interfaceStandard(obj,boundary,neighbour_scheme,neighbour_boundary,type)
            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain

            % Get boundary operators
            [e_neighbour, d_neighbour] = neighbour_scheme.get_boundary_ops(neighbour_boundary);
            [e, d, H_gamma] = obj.get_boundary_ops(boundary);
            Hi = obj.Hi;
            a = obj.a;

            % Get outward unit normal component
            [~, n] = obj.get_boundary_number(boundary);

            Hi = obj.Hi;
            sigma = -n*1i*a/2;
            tau = -n*(1i*a)'/2;

            closure = tau*Hi*d*H_gamma*e' + sigma*Hi*e*H_gamma*d';
            penalty = -tau*Hi*d*H_gamma*e_neighbour' ...
                      -sigma*Hi*e*H_gamma*d_neighbour';

        end

        function [closure, penalty] = interfaceNonConforming(obj,boundary,neighbour_scheme,neighbour_boundary,type)

            % User can request special interpolation operators by specifying type.interpOpSet
            default_field(type, 'interpOpSet', @sbp.InterpOpsOP);
            interpOpSet = type.interpOpSet;

            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            [e_v, d_v] = neighbour_scheme.get_boundary_ops(neighbour_boundary);
            [e_u, d_u, H_gamma] = obj.get_boundary_ops(boundary);
            Hi = obj.Hi;
            a = obj.a;

            % Get outward unit normal component
            [~, n] = obj.get_boundary_number(boundary);


            % Find the number of grid points along the interface
            m_u = size(e_u, 2);
            m_v = size(e_v, 2);

            % Build interpolation operators
            intOps = interpOpSet(m_u, m_v, obj.order, neighbour_scheme.order);
            Iu2v = intOps.Iu2v;
            Iv2u = intOps.Iv2u;

            sigma = -n*1i*a/2;
            tau = -n*(1i*a)'/2;

            closure = tau*Hi*d_u*H_gamma*e_u' + sigma*Hi*e_u*H_gamma*d_u';
            penalty = -tau*Hi*d_u*H_gamma*Iv2u.good*e_v' ...
                      -sigma*Hi*e_u*H_gamma*Iv2u.bad*d_v';

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

        % Returns the boundary ops and sign for the boundary specified by the string boundary.
        % The right boundary is considered the positive boundary
        function [e, d, H_b] = get_boundary_ops(obj, boundary)

            switch boundary
                case 'w'
                    e = obj.e_w;
                    d = obj.d_w;
                    H_b = obj.H_boundary{1};
                case 'e'
                    e = obj.e_e;
                    d = obj.d_e;
                    H_b = obj.H_boundary{1};
                case 's'
                    e = obj.e_s;
                    d = obj.d_s;
                    H_b = obj.H_boundary{2};
                case 'n'
                    e = obj.e_n;
                    d = obj.d_n;
                    H_b = obj.H_boundary{2};
                otherwise
                    error('No such boundary: boundary = %s',boundary);
            end
        end

        function N = size(obj)
            N = prod(obj.m);
        end
    end
end
