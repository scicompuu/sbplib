classdef Hypsyst2d < scheme.Scheme
    properties
        m % Number of points in each direction, possibly a vector
        n %size of system
        h % Grid spacing
        x,y % Grid
        X,Y % Values of x and y for each grid point
        order % Order accuracy for the approximation

        D % non-stabalized scheme operator
        A, B, E %Coefficient matrices

        H % Discrete norm
        % Norms in the x and y directions
        Hxi,Hyi % Kroneckerd norms. 1'*Hx*v corresponds to integration in the x dir.
        I_x,I_y, I_N
        e_w, e_e, e_s, e_n
        params %parameters for the coeficient matrice
    end

    methods
        %Solving Hyperbolic systems on the form u_t=-Au_x-Bu_y-Eu
        function obj = Hypsyst2d(m, lim, order, A, B, E, params)
            default_arg('E', [])
            xlim = lim{1};
            ylim = lim{2};

            if length(m) == 1
                m = [m m];
            end

            obj.A=A;
            obj.B=B;
            obj.E=E;

            m_x = m(1);
            m_y = m(2);
            obj.params = params;

            ops_x = sbp.D2Standard(m_x,xlim,order);
            ops_y = sbp.D2Standard(m_y,ylim,order);

            obj.x = ops_x.x;
            obj.y = ops_y.x;

            obj.X = kr(obj.x,ones(m_y,1));
            obj.Y = kr(ones(m_x,1),obj.y);

            Aevaluated = obj.evaluateCoefficientMatrix(A, obj.X, obj.Y);
            Bevaluated = obj.evaluateCoefficientMatrix(B, obj.X, obj.Y);
            Eevaluated = obj.evaluateCoefficientMatrix(E, obj.X, obj.Y);

            obj.n = length(A(obj.params,0,0));

            I_n = eye(obj.n);I_x = speye(m_x);
            obj.I_x = I_x;
            I_y = speye(m_y);
            obj.I_y = I_y;


            D1_x = kr(I_n, ops_x.D1, I_y);
            obj.Hxi = kr(I_n, ops_x.HI, I_y);
            D1_y = kr(I_n, I_x, ops_y.D1);
            obj.Hyi = kr(I_n, I_x, ops_y.HI);

            obj.e_w = kr(I_n, ops_x.e_l, I_y);
            obj.e_e = kr(I_n, ops_x.e_r, I_y);
            obj.e_s = kr(I_n, I_x, ops_y.e_l);
            obj.e_n = kr(I_n, I_x, ops_y.e_r);

            obj.m = m;
            obj.h = [ops_x.h ops_y.h];
            obj.order = order;

            obj.D = -Aevaluated*D1_x-Bevaluated*D1_y-Eevaluated;

        end

        % Closure functions return the opertors applied to the own doamin to close the boundary
        % Penalty functions return the opertors to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of boundary condition if there are several.
        %       data                is a function returning the data that should be applied at the boundary.
        function [closure, penalty] = boundary_condition(obj,boundary,type,L)
            default_arg('type','char');
            switch type
                case{'c','char'}
                    [closure,penalty] = boundary_condition_char(obj,boundary);
                case{'general'}
                    [closure,penalty] = boundary_condition_general(obj,boundary,L);
                otherwise
                    error('No such boundary condition')
            end
        end

        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary,type)
            error('An interface function does not exist yet');
        end

        function N = size(obj)
            N = obj.m;
        end

        function [ret] = evaluateCoefficientMatrix(obj, mat, X, Y)
            params = obj.params;

            if isa(mat,'function_handle')
                [rows,cols] = size(mat(params,0,0));
                matVec = mat(params,X',Y');
                matVec = sparse(matVec);
                side = max(length(X),length(Y));
            else
                matVec = mat;
                [rows,cols] = size(matVec);
                side = max(length(X),length(Y));
                cols = cols/side;
            end
            ret = cell(rows,cols);

            for ii = 1:rows
                for jj=1:cols
                    ret{ii,jj} = diag(matVec(ii,(jj-1)*side+1:jj*side));
                end
            end
            ret = cell2mat(ret);
        end

        %Characteristic boundary conditions
        function [closure, penalty] = boundary_condition_char(obj,boundary)
            params = obj.params;
            x = obj.x;
            y = obj.y;

            switch boundary
                case {'w','W','west'}
                    e_ = obj.e_w;
                    mat = obj.A;
                    boundPos = 'l';
                    Hi = obj.Hxi;
                    [V,Vi,D,signVec] = obj.matrixDiag(mat,x(1),y);
                    side = max(length(y));
                case {'e','E','east'}
                    e_ = obj.e_e;
                    mat = obj.A;
                    boundPos = 'r';
                    Hi = obj.Hxi;
                    [V,Vi,D,signVec] = obj.matrixDiag(mat,x(end),y);
                    side = max(length(y));
                case {'s','S','south'}
                    e_ = obj.e_s;
                    mat = obj.B;
                    boundPos = 'l';
                    Hi = obj.Hyi;
                    [V,Vi,D,signVec] = obj.matrixDiag(mat,x,y(1));
                    side = max(length(x));
                case {'n','N','north'}
                    e_ = obj.e_n;
                    mat = obj.B;
                    boundPos = 'r';
                    Hi = obj.Hyi;
                    [V,Vi,D,signVec] = obj.matrixDiag(mat,x,y(end));
                    side = max(length(x));
            end
            pos = signVec(1);
            zeroval = signVec(2);
            neg = signVec(3);

            switch boundPos
                case {'l'}
                    tau = sparse(obj.n*side,pos);
                    Vi_plus = Vi(1:pos,:);
                    tau(1:pos,:) = -abs(D(1:pos,1:pos));
                    closure = Hi*e_*V*tau*Vi_plus*e_';
                    penalty = -Hi*e_*V*tau*Vi_plus;
                case {'r'}
                    tau = sparse(obj.n*side,neg);
                    tau((pos+zeroval)+1:obj.n*side,:) = -abs(D((pos+zeroval)+1:obj.n*side,(pos+zeroval)+1:obj.n*side));
                    Vi_minus = Vi((pos+zeroval)+1:obj.n*side,:);
                    closure = Hi*e_*V*tau*Vi_minus*e_';
                    penalty = -Hi*e_*V*tau*Vi_minus;
            end
        end

        % General boundary condition in the form Lu=g(x)
        function [closure,penalty] = boundary_condition_general(obj,boundary,L)
            params = obj.params;
            x = obj.x;
            y = obj.y;

            switch boundary
                case {'w','W','west'}
                    e_ = obj.e_w;
                    mat = obj.A;
                    boundPos = 'l';
                    Hi = obj.Hxi;
                    [V,Vi,D,signVec] = obj.matrixDiag(mat,x(1),y);
                    L = obj.evaluateCoefficientMatrix(L,x(1),y);
                    side = max(length(y));
                case {'e','E','east'}
                    e_ = obj.e_e;
                    mat = obj.A;
                    boundPos = 'r';
                    Hi = obj.Hxi;
                    [V,Vi,D,signVec] = obj.matrixDiag(mat,x(end),y);
                    L = obj.evaluateCoefficientMatrix(L,x(end),y);
                    side = max(length(y));
                case {'s','S','south'}
                    e_ = obj.e_s;
                    mat = obj.B;
                    boundPos = 'l';
                    Hi = obj.Hyi;
                    [V,Vi,D,signVec] = obj.matrixDiag(mat,x,y(1));
                    L = obj.evaluateCoefficientMatrix(L,x,y(1));
                    side = max(length(x));
                case {'n','N','north'}
                    e_ = obj.e_n;
                    mat = obj.B;
                    boundPos = 'r';
                    Hi = obj.Hyi;
                    [V,Vi,D,signVec] = obj.matrixDiag(mat,x,y(end));
                    L = obj.evaluateCoefficientMatrix(L,x,y(end));
                    side = max(length(x));
            end

            pos = signVec(1);
            zeroval = signVec(2);
            neg = signVec(3);

            switch boundPos
                case {'l'}
                    tau = sparse(obj.n*side,pos);
                    Vi_plus = Vi(1:pos,:);
                    Vi_minus = Vi(pos+zeroval+1:obj.n*side,:);
                    V_plus = V(:,1:pos);
                    V_minus = V(:,(pos+zeroval)+1:obj.n*side);

                    tau(1:pos,:) = -abs(D(1:pos,1:pos));
                    R = -inv(L*V_plus)*(L*V_minus);
                    closure = Hi*e_*V*tau*(Vi_plus-R*Vi_minus)*e_';
                    penalty = -Hi*e_*V*tau*inv(L*V_plus)*L;
                case {'r'}
                    tau = sparse(obj.n*side,neg);
                    tau((pos+zeroval)+1:obj.n*side,:) = -abs(D((pos+zeroval)+1:obj.n*side,(pos+zeroval)+1:obj.n*side));
                    Vi_plus = Vi(1:pos,:);
                    Vi_minus = Vi((pos+zeroval)+1:obj.n*side,:);

                    V_plus = V(:,1:pos);
                    V_minus = V(:,(pos+zeroval)+1:obj.n*side);
                    R = -inv(L*V_minus)*(L*V_plus);
                    closure = Hi*e_*V*tau*(Vi_minus-R*Vi_plus)*e_';
                    penalty = -Hi*e_*V*tau*inv(L*V_minus)*L;
            end
        end

        % Function that diagonalizes a symbolic matrix A as A=V*D*Vi
        % D         is a diagonal matrix with the eigenvalues on A on the diagonal sorted by their sign
        %                                    [d+       ]
        %                               D =  [   d0    ]
        %                                    [       d-]
        % signVec   is a vector specifying the number of possitive, zero and negative eigenvalues of D
        function [V,Vi, D,signVec] = matrixDiag(obj,mat,x,y)
            params = obj.params;
            syms xs ys
            [V, D]= eig(mat(params,xs,ys));
            Vi = inv(V);
            xs = x;
            ys = y;

            side = max(length(x),length(y));
            Dret = zeros(obj.n,side*obj.n);
            Vret = zeros(obj.n,side*obj.n);
            Viret = zeros(obj.n,side*obj.n);

            for ii = 1:obj.n
                for jj = 1:obj.n
                    Dret(jj,(ii-1)*side+1:side*ii) = eval(D(jj,ii));
                    Vret(jj,(ii-1)*side+1:side*ii) = eval(V(jj,ii));
                    Viret(jj,(ii-1)*side+1:side*ii) = eval(Vi(jj,ii));
                end
            end

            D = sparse(Dret);
            V = sparse(Vret);
            Vi = sparse(Viret);
            V = obj.evaluateCoefficientMatrix(V,x,y);
            Vi = obj.evaluateCoefficientMatrix(Vi,x,y);
            D = obj.evaluateCoefficientMatrix(D,x,y);
            DD = diag(D);

            poseig = (DD>0);
            zeroeig = (DD==0);
            negeig = (DD<0);

            D = diag([DD(poseig); DD(zeroeig); DD(negeig)]);
            V = [V(:,poseig) V(:,zeroeig) V(:,negeig)];
            Vi = [Vi(poseig,:); Vi(zeroeig,:); Vi(negeig,:)];
            signVec = [sum(poseig),sum(zeroeig),sum(negeig)];
        end

    end
end