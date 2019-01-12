classdef Hypsyst2dCurve < scheme.Scheme
    properties
        m % Number of points in each direction, possibly a vector
        n % size of system
        h % Grid spacing
        X,Y % Values of x and y for each grid point

        J, Ji % Jacobaian and inverse Jacobian
        xi,eta
        Xi,Eta

        A,B
        X_eta, Y_eta
        X_xi,Y_xi
        order % Order accuracy for the approximation

        D % non-stabalized scheme operator
        Ahat, Bhat, E

        H % Discrete norm
        Hxii,Hetai % Kroneckerd norms in xi and eta.
        I_xi,I_eta, I_N, onesN
        e_w, e_e, e_s, e_n
        index_w, index_e,index_s,index_n
        params % Parameters for the coeficient matrice
    end


    methods
        % Solving Hyperbolic systems on the form u_t=-Au_x-Bu_y-Eu
        function obj = Hypsyst2dCurve(m, order, A, B, E, params, ti)
            default_arg('E', [])
            xilim = {0 1};
            etalim = {0 1};

            if length(m) == 1
                m = [m m];
            end
            obj.params = params;
            obj.A=A;
            obj.B=B;

            obj.Ahat=@(params,x,y,x_eta,y_eta)(A(params,x,y).*y_eta-B(params,x,y).*x_eta);
            obj.Bhat=@(params,x,y,x_xi,y_xi)(B(params,x,y).*x_xi-A(params,x,y).*y_xi);
            obj.E=@(params,x,y,~,~)E(params,x,y);

            m_xi = m(1);
            m_eta = m(2);
            m_tot=m_xi*m_eta;

            ops_xi = sbp.D2Standard(m_xi,xilim,order);
            ops_eta = sbp.D2Standard(m_eta,etalim,order);

            obj.xi = ops_xi.x;
            obj.eta = ops_eta.x;

            obj.Xi = kr(obj.xi,ones(m_eta,1));
            obj.Eta = kr(ones(m_xi,1),obj.eta);

            obj.n = length(A(obj.params,0,0));
            obj.onesN=ones(obj.n);

            obj.index_w=1:m_eta;
            obj.index_e=(m_tot-m_e

        metric_termsta+1):m_tot;
            obj.index_s=1:m_eta:(m_tot-m_eta+1);
            obj.index_n=(m_eta):m_eta:m_tot;

            I_n = eye(obj.n);
            I_xi = speye(m_xi);
            obj.I_xi = I_xi;
            I_eta = speye(m_eta);
            obj.I_eta = I_eta;

            D1_xi = kr(I_n, ops_xi.D1, I_eta);
            obj.Hxii = kr(I_n, ops_xi.HI, I_eta);
            D1_eta = kr(I_n, I_xi, ops_eta.D1);
            obj.Hetai = kr(I_n, I_xi, ops_eta.HI);

            obj.e_w = kr(I_n, ops_xi.e_l, I_eta);
            obj.e_e = kr(I_n, ops_xi.e_r, I_eta);
            obj.e_s = kr(I_n, I_xi, ops_eta.e_l);
            obj.e_n = kr(I_n, I_xi,

        metric_termsops_eta.e_r);

            [X,Y] = ti.map(obj.xi,obj.eta);

            [x_xi,x_eta] = gridDerivatives(X,ops_xi.D1,ops_eta.D1);
            [y_xi,y_eta] = gridDerivatives(Y,ops_xi.D1, ops_eta.D1);

            obj.X = reshape(X,m_tot,1);
            obj.Y = reshape(Y,m_tot,1);
            obj.X_xi = reshape(x_xi,m_tot,1);
            obj.Y_xi = reshape(y_xi,m_tot,1);
            obj.X_eta = reshape(x_eta,m_tot,1);
            obj.Y_eta = reshape(y_eta,m_tot,1);

            Ahat_evaluated = obj.evaluateCoefficientMatrix(obj.Ahat, obj.X, obj.Y,obj.X_eta,obj.Y_eta);
            Bhat_evaluated = obj.evaluateCoefficientMatrix(obj.Bhat, obj.X, obj.Y,obj.X_xi,obj.Y_xi);
            E_evaluated = obj.evaluateCoefficientMatrix(obj.E, obj.X, obj.Y,[],[]);

            obj.m = m;
            obj.h = [ops_xi.h ops_eta.h];
            obj.order = order;
            obj.J = obj.X_xi.*obj.Y_eta - obj.X_eta.*obj.Y_xi;
            obj.Ji = kr(I_n,spdiags(1./obj.J,0,m_tot,m_tot));

            obj.D = obj.Ji*(-Ahat_evaluated*D1_xi-Bhat_evaluated*D1_eta)-E_evaluated;

        end

        % Closure functions return the opertors applied to the own doamin to close the boundary
        % Penalty functions return the opertors to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w',General boundary conditions'n','s'.
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

        function [closure, penalty] = interface(obj, boundary, neighbour_scheme, neighbour_boundary, type)
            error('Not implemented');
        end

        function N = size(obj)
            N = obj.m;
        end

        function [ret] = evaluateCoefficientMatrix(obj, mat, X, Y,x_,y_)
            params = obj.params;

            if isa(mat,'function_handle')
                [rows,cols] = size(mat(params,0,0,0,0));
                x_ = kr(obj.onesN,x_);
                y_ = kr(obj.onesN,y_);
                matVec = mat(params,X',Y',x_',y_');
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
                for jj = 1:cols
                    ret{ii,jj} = diag(matVec(ii,(jj-1)*side+1:jj*side));
                end
            end
            ret = cell2mat(ret);
        end

        %Characteristic boundary conditions
        function [closure, penalty] = boundary_condition_char(obj,boundary)
            params = obj.params;
            X = obj.X;
            Y = obj.Y;
            xi = obj.xi;
            eta = obj.eta;
            e_ = obj.getBoundaryOperator('e', boundary);

            switch boundary
                case {'w','W','west'}
                    mat = obj.Ahat;
                    boundPos = 'l';
                    Hi = obj.Hxii;
                    [V,Vi,D,signVec] = obj.matrixDiag(mat,X(obj.index_w),Y(obj.index_w),obj.X_eta(obj.index_w),obj.Y_eta(obj.index_w));
                    side = max(length(eta));
                case {'e','E','east'}
                    mat = obj.Ahat;
                    boundPos = 'r';
                    Hi = obj.Hxii;
                    [V,Vi,D,signVec] = obj.matrixDiag(mat,X(obj.index_e),Y(obj.index_e),obj.X_eta(obj.index_e),obj.Y_eta(obj.index_e));
                    side = max(length(eta));
                case {'s','S','south'}
                    mat = obj.Bhat;
                    boundPos = 'l';
                    Hi = obj.Hetai;
                    [V,Vi,D,signVec] = obj.matrixDiag(mat,X(obj.index_s),Y(obj.index_s),obj.X_xi(obj.index_s),obj.Y_xi(obj.index_s));
                    side = max(length(xi));
                case {'n','N','north'}
                    mat = obj.Bhat;
                    boundPos = 'r';
                    Hi = obj.Hetai;
                    [V,Vi,D,signVec] = obj.matrixDiag(mat,X(obj.index_n),Y(obj.index_n),obj.X_xi(obj.index_n),obj.Y_xi(obj.index_n));
                    side = max(length(xi));
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
            X = obj.X;
            Y = obj.Y;
            xi = obj.xi;
            eta = obj.eta;

            switch boundary
                case {'w','W','west'}
                    e_ = obj.e_w;
                    mat = obj.Ahat;
                    boundPos = 'l';
                    Hi = obj.Hxii;
                    [V,Vi,D,signVec] = obj.matrixDiag(mat,X(obj.index_w),Y(obj.index_w),obj.X_eta(obj.index_w),obj.Y_eta(obj.index_w));

                    Ji_vec = diag(obj.Ji);
                    Ji = diag(Ji_vec(obj.index_w));
                    xi_x = Ji*obj.Y_eta(obj.index_w);
                    xi_y = -Ji*obj.X_eta(obj.index_w);
                    L = obj.evaluateCoefficientMatrix(L,xi_x,xi_y,[],[]);
                    side = max(length(eta));
                case {'e','E','east'}
                    e_ = obj.e_e;
                    mat = obj.Ahat;
                    boundPos = 'r';
                    Hi = obj.Hxii;
                    [V,Vi,D,signVec] = obj.matrixDiag(mat,X(obj.index_e),Y(obj.index_e),obj.X_eta(obj.index_e),obj.Y_eta(obj.index_e));

                    Ji_vec = diag(obj.Ji);
                    Ji = diag(Ji_vec(obj.index_e));
                    xi_x = Ji*obj.Y_eta(obj.index_e);
                    xi_y = -Ji*obj.X_eta(obj.index_e);
                    L = obj.evaluateCoefficientMatrix(L,-xi_x,-xi_y,[],[]);
                    side = max(length(eta));
                case {'s','S','south'}
                    e_ = obj.e_s;
                    mat = obj.Bhat;
                    boundPos = 'l';
                    Hi = obj.Hetai;
                    [V,Vi,D,signVec] = obj.matrixDiag(mat,X(obj.index_s),Y(obj.index_s),obj.X_xi(obj.index_s),obj.Y_xi(obj.index_s));

                    Ji_vec = diag(obj.Ji);
                    Ji = diag(Ji_vec(obj.index_s));
                    eta_x = Ji*obj.Y_xi(obj.index_s);
                    eta_y = -Ji*obj.X_xi(obj.index_s);
                    L = obj.evaluateCoefficientMatrix(L,eta_x,eta_y,[],[]);
                    side = max(length(xi));
                case {'n','N','north'}
                    e_ = obj.e_n;
                    mat = obj.Bhat;
                    boundPos = 'r';
                    Hi = obj.Hetai;
                    [V,Vi,D,signVec] = obj.matrixDiag(mat,X(obj.index_n),Y(obj.index_n),obj.X_xi(obj.index_n),obj.Y_xi(obj.index_n));

                    Ji_vec = diag(obj.Ji);
                    Ji = diag(Ji_vec(obj.index_n));
                    eta_x = Ji*obj.Y_xi(obj.index_n);
                    eta_y = -Ji*obj.X_xi(obj.index_n);
                    L = obj.evaluateCoefficientMatrix(L,-eta_x,-eta_y,[],[]);
                    side = max(length(xi));
            end

            pos = signVec(1);
            zeroval = signVec(2);
            neg = signVec(3);

            switch boundPos
                case {'l'}
                    tau = sparse(obj.n*side,pos);
                    Vi_plus = Vi(1:pos,:);
                    Vi_minus = Vi(pos+1:obj.n*side,:);
                    V_plus = V(:,1:pos);
                    V_minus = V(:,(pos)+1:obj.n*side);

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
        function [V,Vi, D,signVec] = matrixDiag(obj,mat,x,y,x_,y_)
            params = obj.params;
            syms xs ys
            if(sum(abs(x_)) ~= 0)
                syms xs_
            else
                xs_ = 0;
            end

            if(sum(abs(y_))~= 0)
                syms ys_;
            else
                ys_ = 0;
            end

            [V, D] = eig(mat(params,xs,ys,xs_,ys_));
            Vi = inv(V);
            syms xs ys xs_ ys_

            xs = x;
            ys = y;
            xs_ = x_;
            ys_ = y_;

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
            V = obj.evaluateCoefficientMatrix(V,x,y,x_,y_);
            D = obj.evaluateCoefficientMatrix(D,x,y,x_,y_);
            Vi = obj.evaluateCoefficientMatrix(Vi,x,y,x_,y_);
            DD = diag(D);

            poseig = (DD>0);
            zeroeig = (DD==0);
            negeig = (DD<0);

            D = diag([DD(poseig); DD(zeroeig); DD(negeig)]);
            V = [V(:,poseig) V(:,zeroeig) V(:,negeig)];
            Vi = [Vi(poseig,:); Vi(zeroeig,:); Vi(negeig,:)];
            signVec = [sum(poseig),sum(zeroeig),sum(negeig)];
        end

        % Returns the boundary operator op for the boundary specified by the string boundary.
        % op        -- string or a cell array of strings
        % boundary  -- string
        function varargout = getBoundaryOperator(obj, op, boundary)

            if ~iscell(op)
                op = {op};
            end

            for i = 1:numel(op)
                switch op{i}
                case 'e'
                    switch boundary
                    case 'w'
                        e = obj.e_w;
                    case 'e'
                        e = obj.e_e;
                    case 's'
                        e = obj.e_s;
                    case 'n'
                        e = obj.e_n;
                    otherwise
                        error('No such boundary: boundary = %s',boundary);
                    end
                    varargout{i} = e;
                end
            end
        end

        % Returns square boundary quadrature matrix, of dimension
        % corresponding to the number of boundary points
        %
        % boundary -- string
        function H_b = getBoundaryQuadrature(obj, boundary)

            e = obj.getBoundaryOperator('e', boundary);

            switch boundary
                case 'w'
                    H_b = inv(e'*obj.Hetai*e);
                case 'e'
                    H_b = inv(e'*obj.Hetai*e);
                case 's'
                    H_b = inv(e'*obj.Hxii*e);
                case 'n'
                    H_b = inv(e'*obj.Hxii*e);
                otherwise
                    error('No such boundary: boundary = %s',boundary);
            end
        end


    end
end