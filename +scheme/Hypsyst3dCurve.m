classdef Hypsyst3dCurve < scheme.Scheme
    properties
        m % Number of points in each direction, possibly a vector
        n %size of system
        h % Grid spacing
        X, Y, Z% Values of x and y for each grid point
        Yx, Zx, Xy, Zy, Xz, Yz %Grid values for boundary surfaces
        
        xi,eta,zeta
        Xi, Eta, Zeta
        
        Eta_xi, Zeta_xi, Xi_eta, Zeta_eta, Xi_zeta, Eta_zeta
        
        X_xi, X_eta, X_zeta,Y_xi,Y_eta,Y_zeta,Z_xi,Z_eta,Z_zeta
        Aev
        
        metric_terms
        
        order % Order accuracy for the approximation
        
        D % non-stabalized scheme operator
        Aevaluated, Bevaluated, Cevaluated, Eevaluated
        Ahat, Bhat, Chat, E
        A,B,C
        
        J, Ji
        
        H % Discrete norm
        % Norms in the x, y and z directions
        Hxii,Hetai,Hzetai, Hzi % Kroneckerd norms. 1'*Hx*v corresponds to integration in the x dir.
        Hxi,Heta,Hzeta
        I_xi,I_eta,I_zeta, I_N,onesN
        e_w, e_e, e_s, e_n, e_b, e_t
        index_w, index_e,index_s,index_n, index_b, index_t
        params %parameters for the coeficient matrice
    end
    
    
    methods
        function obj = Hypsyst3dCurve(m, order, A, B,C, E, params,ti,operator)
            xilim ={0 1};
            etalim = {0 1};
            zetalim = {0 1};
            
            if length(m) == 1
                m = [m m m];
            end
            m_xi = m(1);
            m_eta = m(2);
            m_zeta = m(3);
            m_tot = m_xi*m_eta*m_zeta;
            obj.params = params;
            obj.n = length(A(obj,0,0,0));
            
            obj.m = m;
            obj.order = order;
            obj.onesN = ones(obj.n);
            
            switch operator
                case 'upwind'
                    ops_xi = sbp.D1Upwind(m_xi,xilim,order);
                    ops_eta = sbp.D1Upwind(m_eta,etalim,order);
                    ops_zeta = sbp.D1Upwind(m_zeta,zetalim,order);
                case 'standard'
                    ops_xi = sbp.D2Standard(m_xi,xilim,order);
                    ops_eta = sbp.D2Standard(m_eta,etalim,order);
                    ops_zeta = sbp.D2Standard(m_zeta,zetalim,order); 
                otherwise
                    error('Operator not available')
            end
            
            obj.xi = ops_xi.x;
            obj.eta = ops_eta.x;
            obj.zeta = ops_zeta.x;
            
            obj.Xi = kr(obj.xi,ones(m_eta,1),ones(m_zeta,1));%% Que pasa?
            obj.Eta = kr(ones(m_xi,1),obj.eta,ones(m_zeta,1));
            obj.Zeta = kr(ones(m_xi,1),ones(m_eta,1),obj.zeta);
            
            obj.Eta_xi = kr(obj.eta,ones(m_xi,1));
            obj.Zeta_xi = kr(ones(m_eta,1),obj.zeta);
            
            obj.Xi_eta = kr(obj.xi,ones(m_zeta,1));
            obj.Zeta_eta = kr(ones(m_xi,1),obj.zeta);
            
            obj.Xi_zeta = kr(obj.xi,ones(m_eta,1));
            obj.Eta_zeta = kr(ones(m_zeta,1),obj.eta);
            
            [X,Y,Z] = ti.map(obj.Xi,obj.Eta,obj.Zeta);
            obj.X = X;
            obj.Y = Y;
            obj.Z = Z;
            
            I_n = eye(obj.n);
            I_xi = speye(m_xi);
            obj.I_xi = I_xi;
            I_eta = speye(m_eta);
            obj.I_eta = I_eta;
            I_zeta = speye(m_zeta);
            obj.I_zeta = I_zeta;     
            
            I_N=kr(I_n,I_xi,I_eta,I_zeta);
            
            O_xi = ones(m_xi,1);
            O_eta = ones(m_eta,1);
            O_zeta = ones(m_zeta,1);
            
            obj.Hxii = kr(I_n, ops_xi.HI, I_eta,I_zeta);
            obj.Hetai = kr(I_n, I_xi, ops_eta.HI,I_zeta);           
            obj.Hzetai = kr(I_n, I_xi,I_eta, ops_zeta.HI);
            obj.Hxi = ops_xi.H;
            obj.Heta = ops_eta.H;           
            obj.Hzeta = ops_zeta.H;
            obj.h = [ops_xi.h ops_eta.h ops_zeta.h];
            
            switch operator
                case 'upwind'
                D1_xi = kr((ops_xi.Dp+ops_xi.Dm)/2, I_eta,I_zeta);
                D1_eta = kr(I_xi, (ops_eta.Dp+ops_eta.Dm)/2,I_zeta);
                D1_zeta = kr(I_xi, I_eta,(ops_zeta.Dp+ops_zeta.Dm)/2);
                otherwise
                D1_xi = kr(ops_xi.D1, I_eta,I_zeta);
                D1_eta = kr(I_xi, ops_eta.D1,I_zeta);
                D1_zeta = kr(I_xi, I_eta,ops_zeta.D1);                   
            end
            
            obj.e_w = kr(I_n, ops_xi.e_l, I_eta,I_zeta);
            obj.e_e = kr(I_n, ops_xi.e_r, I_eta,I_zeta);
            obj.e_s = kr(I_n, I_xi, ops_eta.e_l,I_zeta);
            obj.e_n = kr(I_n, I_xi, ops_eta.e_r,I_zeta);
            obj.e_b = kr(I_n, I_xi, I_eta, ops_zeta.e_l);
            obj.e_t = kr(I_n, I_xi, I_eta, ops_zeta.e_r);
            
            obj.A = A;
            obj.B = B;
            obj.C = C;
            
            obj.X_xi = D1_xi*X;
            obj.X_eta = D1_eta*X;
            obj.X_zeta = D1_zeta*X;
            obj.Y_xi = D1_xi*Y;
            obj.Y_eta = D1_eta*Y;
            obj.Y_zeta = D1_zeta*Y;
            obj.Z_xi = D1_xi*Z;
            obj.Z_eta = D1_eta*Z;
            obj.Z_zeta = D1_zeta*Z;
            
            D1_xi = kr(I_n,D1_xi);
            D1_eta = kr(I_n,D1_eta);
            D1_zeta = kr(I_n,D1_zeta);
            
            obj.index_w = (kr(ops_xi.e_l, O_eta,O_zeta)==1);
            obj.index_e = (kr(ops_xi.e_r, O_eta,O_zeta)==1);
            obj.index_s = (kr(O_xi, ops_eta.e_l,O_zeta)==1);
            obj.index_n = (kr(O_xi, ops_eta.e_r,O_zeta)==1);
            obj.index_b = (kr(O_xi, O_eta, ops_zeta.e_l)==1);
            obj.index_t = (kr(O_xi, O_eta, ops_zeta.e_r)==1);
            
            
            obj.Ahat = @transform_coefficient_matrix;
            obj.Bhat = @transform_coefficient_matrix;
            obj.Chat = @transform_coefficient_matrix;
            obj.E = @(obj,x,y,z,~,~,~,~,~,~)E(obj,x,y,z);
            
            obj.Aevaluated = obj.evaluateCoefficientMatrix(obj.Ahat,obj.X, obj.Y,obj.Z, obj.X_eta,obj.X_zeta,obj.Y_eta,obj.Y_zeta,obj.Z_eta,obj.Z_zeta);
            obj.Bevaluated = obj.evaluateCoefficientMatrix(obj.Bhat,obj.X, obj.Y,obj.Z, obj.X_zeta,obj.X_xi,obj.Y_zeta,obj.Y_xi,obj.Z_zeta,obj.Z_xi);
            obj.Cevaluated = obj.evaluateCoefficientMatrix(obj.Chat,obj.X,obj.Y,obj.Z, obj.X_xi,obj.X_eta,obj.Y_xi,obj.Y_eta,obj.Z_xi,obj.Z_eta);
            obj.Eevaluated = obj.evaluateCoefficientMatrix(obj.E, obj.X, obj.Y,obj.Z,[],[],[],[],[],[]);
            
            obj.J = obj.X_xi.*obj.Y_eta.*obj.Z_zeta...
                 +obj.X_zeta.*obj.Y_xi.*obj.Z_eta...
                 +obj.X_eta.*obj.Y_zeta.*obj.Z_xi...
                 -obj.X_xi.*obj.Y_zeta.*obj.Z_eta...
                 -obj.X_eta.*obj.Y_xi.*obj.Z_zeta...
                 -obj.X_zeta.*obj.Y_eta.*obj.Z_xi;
            
            obj.Ji = kr(I_n,spdiags(1./obj.J,0,m_tot,m_tot));
            
            
            switch operator
                case 'upwind'
                    alphaA = max(eig(obj.Ahat(obj,obj.X(end), obj.Y(end),obj.Z(end), obj.X_eta(end),obj.X_zeta(end),obj.Y_eta(end),obj.Y_zeta(end),obj.Z_eta(end),obj.Z_zeta(end))));
                    alphaB = max(eig(obj.Bhat(obj,obj.X(end), obj.Y(end),obj.Z(end), obj.X_zeta(end),obj.X_xi(end),obj.Y_zeta(end),obj.Y_xi(end),obj.Z_zeta(end),obj.Z_xi(end))));
                    alphaC = max(eig(obj.Chat(obj,obj.X(end), obj.Y(end),obj.Z(end), obj.X_xi(end),obj.X_eta(end),obj.Y_xi(end),obj.Y_eta(end),obj.Z_xi(end),obj.Z_eta(end))));
                    
                    Ap = (obj.Aevaluated+alphaA*I_N)/2;
                    Dmxi = kr(I_n, ops_xi.Dm, I_eta,I_zeta);
                    diffSum=-Ap*Dmxi;
                    clear Ap Dmxi

                    Am = (obj.Aevaluated-alphaA*I_N)/2;
                    clear obj.Aevaluated
                    Dpxi = kr(I_n, ops_xi.Dp, I_eta,I_zeta);
                    temp=Am*Dpxi;
                    diffSum=diffSum-temp;
                    clear Am Dpxi

                    Bp = (obj.Bevaluated+alphaB*I_N)/2;
                    Dmeta = kr(I_n, I_xi, ops_eta.Dm,I_zeta);
                    temp=Bp*Dmeta;
                    diffSum=diffSum-temp;
                    clear Bp Dmeta

                    Bm = (obj.Bevaluated-alphaB*I_N)/2;
                    clear obj.Bevaluated
                    Dpeta = kr(I_n, I_xi, ops_eta.Dp,I_zeta);
                    temp=Bm*Dpeta;
                    diffSum=diffSum-temp;


                    Cp = (obj.Cevaluated+alphaC*I_N)/2;
                    Dmzeta = kr(I_n, I_xi, I_eta,ops_zeta.Dm);
                    temp=Cp*Dmzeta;
                    diffSum=diffSum-temp;
                    clear Cp Dmzeta

                    Cm = (obj.Cevaluated-alphaC*I_N)/2;      
                    clear obj.Cevaluated     
                    Dpzeta = kr(I_n, I_xi, I_eta,ops_zeta.Dp);
                    temp=Cm*Dpzeta;
                    diffSum=diffSum-temp;
                  
                    obj.D = obj.Ji*diffSum-obj.Eevaluated;
                otherwise
                    obj.D = obj.Ji*(-obj.Aevaluated*D1_xi-obj.Bevaluated*D1_eta -obj.Cevaluated*D1_zeta)-obj.Eevaluated;
            end
            
            
        end
        
        function [ret] = transform_coefficient_matrix(obj,x,y,z,x_1,x_2,y_1,y_2,z_1,z_2)
            ret = obj.A(obj,x,y,z).*(y_1.*z_2-z_1.*y_2);
            ret = ret+obj.B(obj,x,y,z).*(x_2.*z_1-x_1.*z_2);
            ret = ret+obj.C(obj,x,y,z).*(x_1.*y_2-x_2.*y_1);
        end
        
        
        % Closure functions return the opertors applied to the own doamin to close the boundary
        % Penalty functions return the opertors to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of boundary condition if there are several.
        %       data                is a function returning the data that should be applied at the boundary.
        function [closure, penalty] = boundary_condition(obj,boundary,type,L)
            default_arg('type','char');
            BM = boundary_matrices(obj,boundary);
            
            switch type
                case{'c','char'}
                    [closure,penalty] = boundary_condition_char(obj,BM);
                case{'general'}
                    [closure,penalty] = boundary_condition_general(obj,BM,boundary,L);
                otherwise
                    error('No such boundary condition')
            end
        end
        
        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary)
            error('An interface function does not exist yet');
        end
        
        function N = size(obj)
            N = obj.m;
        end
        
        function [ret] = evaluateCoefficientMatrix(obj,mat, X, Y, Z , x_1 , x_2 , y_1 , y_2 , z_1 , z_2)
            params = obj.params;
            side = max(length(X),length(Y));
            if isa(mat,'function_handle')
                [rows,cols] = size(mat(obj,0,0,0,0,0,0,0,0,0));
                x_1 = kr(obj.onesN,x_1);
                x_2 = kr(obj.onesN,x_2);
                y_1 = kr(obj.onesN,y_1);
                y_2 = kr(obj.onesN,y_2);
                z_1 = kr(obj.onesN,z_1);
                z_2 = kr(obj.onesN,z_2);
                matVec = mat(obj,X',Y',Z',x_1',x_2',y_1',y_2',z_1',z_2');
                matVec = sparse(matVec);
            else
                matVec = mat;
                [rows,cols] = size(matVec);
                side = max(length(X),length(Y));
                cols = cols/side;
            end
            ret = cell(rows,cols);
            
            
            for ii=1:rows
                for jj=1:cols
                    ret{ii,jj} = diag(matVec(ii,(jj-1)*side+1:jj*side));
                end
            end
            
            ret = cell2mat(ret);
        end
        
        
        function [BM] = boundary_matrices(obj,boundary)
            params = obj.params;
            BM.boundary = boundary;
            switch boundary
                case {'w','W','west'}
                    BM.e_ = obj.e_w;
                    mat = obj.Ahat;
                    BM.boundpos = 'l';
                    BM.Hi = obj.Hxii;
                    BM.index = obj.index_w;
                    BM.x_1 = obj.X_eta(BM.index);
                    BM.x_2 = obj.X_zeta(BM.index);
                    BM.y_1 = obj.Y_eta(BM.index);
                    BM.y_2 = obj.Y_zeta(BM.index);
                    BM.z_1 = obj.Z_eta(BM.index);
                    BM.z_2 = obj.Z_zeta(BM.index);
                case {'e','E','east'}
                    BM.e_ = obj.e_e;
                    mat = obj.Ahat;
                    BM.boundpos = 'r';
                    BM.Hi = obj.Hxii;
                    BM.index = obj.index_e;
                    BM.x_1 = obj.X_eta(BM.index);
                    BM.x_2 = obj.X_zeta(BM.index);
                    BM.y_1 = obj.Y_eta(BM.index);
                    BM.y_2 = obj.Y_zeta(BM.index);
                    BM.z_1 = obj.Z_eta(BM.index);
                    BM.z_2 = obj.Z_zeta(BM.index);
                case {'s','S','south'}
                    BM.e_ = obj.e_s;
                    mat = obj.Bhat;
                    BM.boundpos = 'l';
                    BM.Hi = obj.Hetai;
                    BM.index = obj.index_s;
                    BM.x_1 = obj.X_zeta(BM.index);
                    BM.x_2 = obj.X_xi(BM.index);
                    BM.y_1 = obj.Y_zeta(BM.index);
                    BM.y_2 = obj.Y_xi(BM.index);
                    BM.z_1 = obj.Z_zeta(BM.index);
                    BM.z_2 = obj.Z_xi(BM.index);
                case {'n','N','north'}
                    BM.e_ = obj.e_n;
                    mat = obj.Bhat;
                    BM.boundpos = 'r';
                    BM.Hi = obj.Hetai;
                    BM.index = obj.index_n;
                    BM.x_1 = obj.X_zeta(BM.index);
                    BM.x_2 = obj.X_xi(BM.index);
                    BM.y_1 = obj.Y_zeta(BM.index);
                    BM.y_2 = obj.Y_xi(BM.index);
                    BM.z_1 = obj.Z_zeta(BM.index);
                    BM.z_2 = obj.Z_xi(BM.index);
                case{'b','B','Bottom'}
                    BM.e_ = obj.e_b;
                    mat = obj.Chat;
                    BM.boundpos = 'l';
                    BM.Hi = obj.Hzetai;
                    BM.index = obj.index_b;
                    BM.x_1 = obj.X_xi(BM.index);
                    BM.x_2 = obj.X_eta(BM.index);
                    BM.y_1 = obj.Y_xi(BM.index);
                    BM.y_2 = obj.Y_eta(BM.index);
                    BM.z_1 = obj.Z_xi(BM.index);
                    BM.z_2 = obj.Z_eta(BM.index);
                case{'t','T','Top'}
                    BM.e_ = obj.e_t;
                    mat = obj.Chat;
                    BM.boundpos = 'r';
                    BM.Hi = obj.Hzetai;
                    BM.index = obj.index_t;
                    BM.x_1 = obj.X_xi(BM.index);
                    BM.x_2 = obj.X_eta(BM.index);
                    BM.y_1 = obj.Y_xi(BM.index);
                    BM.y_2 = obj.Y_eta(BM.index);
                    BM.z_1 = obj.Z_xi(BM.index);
                    BM.z_2 = obj.Z_eta(BM.index);
            end
            [BM.V,BM.Vi,BM.D,signVec] = obj.matrixDiag(mat,obj.X(BM.index),obj.Y(BM.index),obj.Z(BM.index),...
                BM.x_1,BM.x_2,BM.y_1,BM.y_2,BM.z_1,BM.z_2);
            BM.side = sum(BM.index);
            BM.pos = signVec(1); BM.zeroval=signVec(2); BM.neg=signVec(3);
        end
        
        
        function [closure, penalty]=boundary_condition_char(obj,BM)
            side = BM.side;
            pos = BM.pos;
            neg = BM.neg;
            zeroval = BM.zeroval;
            V = BM.V;
            Vi = BM.Vi;
            Hi = BM.Hi;
            D = BM.D;
            e_ = BM.e_;
            
            switch BM.boundpos
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
        
        
        function [closure,penalty]=boundary_condition_general(obj,BM,boundary,L)
            side = BM.side;
            pos = BM.pos;
            neg = BM.neg;
            zeroval = BM.zeroval;
            V = BM.V;
            Vi = BM.Vi;
            Hi = BM.Hi;
            D = BM.D;
            e_ = BM.e_;
            index = BM.index;
            
            switch BM.boundary
                case{'b','B','bottom'}
                    Ji_vec = diag(obj.Ji);
                    Ji = diag(Ji_vec(index));
                    Zeta_x = Ji*(obj.Y_xi(index).*obj.Z_eta(index)-obj.Z_xi(index).*obj.Y_eta(index));
                    Zeta_y = Ji*(obj.X_eta(index).*obj.Z_xi(index)-obj.X_xi(index).*obj.Z_eta(index));
                    Zeta_z = Ji*(obj.X_xi(index).*obj.Y_eta(index)-obj.Y_xi(index).*obj.X_eta(index));
                    
                    L = obj.evaluateCoefficientMatrix(L,Zeta_x,Zeta_y,Zeta_z,[],[],[],[],[],[]);
            end
            
            switch BM.boundpos
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
        
        
        function [V,Vi, D,signVec] = matrixDiag(obj,mat,x,y,z,x_1,x_2,y_1,y_2,z_1,z_2)
            params = obj.params;
            eps = 10^(-10);
            if(sum(abs(x_1))>eps)
                syms x_1s
            else
                x_1s = 0;
            end
            
            if(sum(abs(x_2))>eps)
                syms x_2s;
            else
                x_2s = 0;
            end
            
            
            if(sum(abs(y_1))>eps)
                syms y_1s
            else
                y_1s = 0;
            end
            
            if(sum(abs(y_2))>eps)
                syms y_2s;
            else
                y_2s = 0;
            end
            
            if(sum(abs(z_1))>eps)
                syms z_1s
            else
                z_1s = 0;
            end
            
            if(sum(abs(z_2))>eps)
                syms z_2s;
            else
                z_2s = 0;
            end
            
            syms xs ys zs
            [V, D] = eig(mat(obj,xs,ys,zs,x_1s,x_2s,y_1s,y_2s,z_1s,z_2s));
            Vi = inv(V);
            %    syms x_1s x_2s y_1s y_2s z_1s z_2s
            xs = x;
            ys = y;
            zs = z;
            x_1s = x_1;
            x_2s = x_2;
            y_1s = y_1;
            y_2s = y_2;
            z_1s = z_1;
            z_2s = z_2;
            
            side = max(length(x),length(y));
            Dret = zeros(obj.n,side*obj.n);
            Vret = zeros(obj.n,side*obj.n);
            Viret = zeros(obj.n,side*obj.n);
            
            for ii=1:obj.n
                for jj=1:obj.n
                    Dret(jj,(ii-1)*side+1:side*ii) = eval(D(jj,ii));
                    Vret(jj,(ii-1)*side+1:side*ii) = eval(V(jj,ii));
                    Viret(jj,(ii-1)*side+1:side*ii) = eval(Vi(jj,ii));
                end
            end
            
            D = sparse(Dret);
            V = sparse(Vret);
            Vi = sparse(Viret);
            V = obj.evaluateCoefficientMatrix(V,x,y,z,x_1,x_2,y_1,y_2,z_1,z_2);
            D = obj.evaluateCoefficientMatrix(D,x,y,z,x_1,x_2,y_1,y_2,z_1,z_2);
            Vi = obj.evaluateCoefficientMatrix(Vi,x,y,z,x_1,x_2,y_1,y_2,z_1,z_2);
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
