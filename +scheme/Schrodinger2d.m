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

        interpolation_type % MC or AWW

    end

    methods

        function obj = Schrodinger2d(g ,order, a, opSet, interpolation_type)
            default_arg('interpolation_type','AWW');
            default_arg('opSet',{@sbp.D2Variable, @sbp.D2Variable});
            default_arg('a',1);
            dim = 2;

            assert(isa(g, 'grid.Cartesian'))

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
            obj.interpolation_type = interpolation_type;

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
            a = obj.a;

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

        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary)
            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            % Get neighbour boundary operator

            [coord_nei, n_nei] = get_boundary_number(obj, neighbour_boundary);
            [coord, n] = get_boundary_number(obj, boundary);
            switch n_nei
            case 1
                % North or east boundary
                e_neighbour = neighbour_scheme.e_r;
                d_neighbour = neighbour_scheme.d1_r;
            case -1
                % South or west boundary
                e_neighbour = neighbour_scheme.e_l;
                d_neighbour = neighbour_scheme.d1_l;
            end

            e_neighbour = e_neighbour{coord_nei};
            d_neighbour = d_neighbour{coord_nei};
            H_gamma = obj.H_boundary{coord};
            Hi = obj.Hi;
            a = obj.a;

            switch coord_nei
            case 1
                m_neighbour = neighbour_scheme.m(2);
            case 2
                m_neighbour = neighbour_scheme.m(1);
            end

            switch coord
            case 1
                m = obj.m(2);
            case 2
                m = obj.m(1);
            end

           switch n
            case 1
                % North or east boundary
                e = obj.e_r;
                d = obj.d1_r;
            case -1
                % South or west boundary
                e = obj.e_l;
                d = obj.d1_l;
            end
            e = e{coord};
            d = d{coord}; 

            Hi = obj.Hi;
            sigma = -n*1i*a/2;
            tau = -n*(1i*a)'/2;

            grid_ratio = m/m_neighbour;
             if grid_ratio ~= 1

                [ms, index] = sort([m, m_neighbour]);
                orders = [obj.order, neighbour_scheme.order];
                orders = orders(index);

                switch obj.interpolation_type
                case 'MC'
                    interpOpSet = sbp.InterpMC(ms(1),ms(2),orders(1),orders(2));
                    if grid_ratio < 1
                        I_neighbour2local_e = interpOpSet.IF2C;
                        I_neighbour2local_d = interpOpSet.IF2C;
                        I_local2neighbour_e = interpOpSet.IC2F;
                        I_local2neighbour_d = interpOpSet.IC2F;
                    elseif grid_ratio > 1
                        I_neighbour2local_e = interpOpSet.IC2F;
                        I_neighbour2local_d = interpOpSet.IC2F;
                        I_local2neighbour_e = interpOpSet.IF2C;
                        I_local2neighbour_d = interpOpSet.IF2C;
                    end
                case 'AWW'
                    %String 'C2F' indicates that ICF2 is more accurate.
                    interpOpSetF2C = sbp.InterpAWW(ms(1),ms(2),orders(1),orders(2),'F2C');
                    interpOpSetC2F = sbp.InterpAWW(ms(1),ms(2),orders(1),orders(2),'C2F'); 
                    if grid_ratio < 1 
                        % Local is coarser than neighbour
                        I_neighbour2local_e = interpOpSetF2C.IF2C;
                        I_neighbour2local_d = interpOpSetC2F.IF2C;
                        I_local2neighbour_e = interpOpSetC2F.IC2F;
                        I_local2neighbour_d = interpOpSetF2C.IC2F;
                    elseif grid_ratio > 1
                        % Local is finer than neighbour 
                        I_neighbour2local_e = interpOpSetC2F.IC2F;
                        I_neighbour2local_d = interpOpSetF2C.IC2F;
                        I_local2neighbour_e = interpOpSetF2C.IF2C;
                        I_local2neighbour_d = interpOpSetC2F.IF2C;
                    end
                otherwise
                    error(['Interpolation type ' obj.interpolation_type ...
                         ' is not available.' ]);
                end

             else 
                % No interpolation required
                I_neighbour2local_e = speye(m,m);
                I_neighbour2local_d = speye(m,m);
                I_local2neighbour_e = speye(m,m);
                I_local2neighbour_d = speye(m,m);
            end

            closure = tau*Hi*d*H_gamma*e' + sigma*Hi*e*H_gamma*d';
            penalty = -tau*Hi*d*H_gamma*I_neighbour2local_e*e_neighbour' ...
                      -sigma*Hi*e*H_gamma*I_neighbour2local_d*d_neighbour'; 
             
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
                    error(['No such operator: operator = ' op]);
            end

        end

        function N = size(obj)
            N = prod(obj.m);
        end
    end
end
