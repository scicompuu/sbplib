classdef Hypsyst3d < scheme.Scheme
    properties
        m % Number of points in each direction, possibly a vector
        n %size of system
        h % Grid spacing
        x, y, z % Grid
        X, Y, Z% Values of x and y for each grid point
        Yx, Zx, Xy, Zy, Xz, Yz %Grid values for boundary surfaces
        order % Order accuracy for the approximation
        
        D % non-stabalized scheme operator
        A, B, C, E
        Aevaluated,Bevaluated,Cevaluated, Eevaluated
        
        H % Discrete norm
        % Norms in the x, y and z directions
        Hx, Hy, Hz
        Hxi,Hyi, Hzi % Kroneckerd norms. 1'*Hx*v corresponds to integration in the x dir.
        I_x,I_y, I_z, I_N
        e_w, e_e, e_s, e_n, e_b, e_t
        params %parameters for the coeficient matrice
    end
    
    
    methods
        function obj = Hypsyst3d(m, lim, order, A, B,C, E, params,operator)
            default_arg('E', [])
            default_arg('operatpr',[])
            xlim = lim{1};
            ylim = lim{2};
            zlim = lim{3};
            
            if length(m) == 1
                m = [m m m];
            end
            
            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.E = E;
            m_x = m(1);
            m_y = m(2);
            m_z=m(3);
            obj.params = params;
            
            switch operator
                case 'upwind'
                    ops_x = sbp.D1Upwind(m_x,xlim,order);
                    ops_y = sbp.D1Upwind(m_y,ylim,order);
                    ops_z = sbp.D1Upwind(m_z,zlim,order);
                otherwise
                    ops_x = sbp.D2Standard(m_x,xlim,order);
                    ops_y = sbp.D2Standard(m_y,ylim,order);
                    ops_z = sbp.D2Standard(m_z,zlim,order);
            end
            
            obj.x = ops_x.x;
            obj.y = ops_y.x;
            obj.z = ops_z.x;
            
            obj.X = kr(obj.x,ones(m_y,1),ones(m_z,1));%% Que pasa?
            obj.Y = kr(ones(m_x,1),obj.y,ones(m_z,1));
            obj.Z = kr(ones(m_x,1),ones(m_y,1),obj.z);
            
            obj.Yx = kr(obj.y,ones(m_z,1));
            obj.Zx = kr(ones(m_y,1),obj.z);
            
            obj.Xy = kr(obj.x,ones(m_z,1));
            obj.Zy = kr(ones(m_x,1),obj.z);
            
            obj.Xz = kr(obj.x,ones(m_y,1));
            obj.Yz = kr(ones(m_z,1),obj.y);
            
            obj.Aevaluated = obj.evaluateCoefficientMatrix(A, obj.X, obj.Y,obj.Z);
            obj.Bevaluated = obj.evaluateCoefficientMatrix(B, obj.X, obj.Y,obj.Z);
            obj.Cevaluated = obj.evaluateCoefficientMatrix(C, obj.X, obj.Y,obj.Z);
            obj.Eevaluated = obj.evaluateCoefficientMatrix(E, obj.X, obj.Y,obj.Z);
            
            obj.n = length(A(obj.params,0,0,0));
            
            I_n = eye(obj.n);
            I_x = speye(m_x);
            obj.I_x = I_x;
            I_y = speye(m_y);
            obj.I_y = I_y;
            I_z = speye(m_z);
            obj.I_z = I_z;
            I_N=kr(I_n,I_x,I_y,I_z);
            
            obj.Hxi = kr(I_n, ops_x.HI, I_y,I_z);
            obj.Hx = ops_x.H;
            obj.Hyi = kr(I_n, I_x, ops_y.HI,I_z);
            obj.Hy = ops_y.H
            obj.Hzi = kr(I_n, I_x,I_y, ops_z.HI);
            obj.Hz = ops_z.H;
            
            obj.e_w = kr(I_n, ops_x.e_l, I_y,I_z);
            obj.e_e = kr(I_n, ops_x.e_r, I_y,I_z);
            obj.e_s = kr(I_n, I_x, ops_y.e_l,I_z);
            obj.e_n = kr(I_n, I_x, ops_y.e_r,I_z);
            obj.e_b = kr(I_n, I_x, I_y, ops_z.e_l);
            obj.e_t = kr(I_n, I_x, I_y, ops_z.e_r);
            
            obj.m = m;
            obj.h = [ops_x.h ops_y.h ops_x.h];
            obj.order = order;
            
            switch operator
                case 'upwind'
                    alphaA = max(eig(A(params,obj.x(end),obj.y(end),obj.z(end))));
                    alphaB = max(eig(B(params,obj.x(end),obj.y(end),obj.z(end))));
                    alphaC = max(eig(C(params,obj.x(end),obj.y(end),obj.z(end))));
                    
                    Ap = (obj.Aevaluated+alphaA*I_N)/2;
                    Am = (obj.Aevaluated-alphaA*I_N)/2;
                    Bp = (obj.Bevaluated+alphaB*I_N)/2;
                    Bm = (obj.Bevaluated-alphaB*I_N)/2;
                    Cp = (obj.Cevaluated+alphaC*I_N)/2;
                    Cm = (obj.Cevaluated-alphaC*I_N)/2;
                    
                    Dpx = kr(I_n, ops_x.Dp, I_y,I_z);
                    Dmx = kr(I_n, ops_x.Dm, I_y,I_z);
                    Dpy = kr(I_n, I_x, ops_y.Dp,I_z);
                    Dmy = kr(I_n, I_x, ops_y.Dm,I_z);
                    Dpz = kr(I_n, I_x, I_y,ops_z.Dp);
                    Dmz = kr(I_n, I_x, I_y,ops_z.Dm);
                    
                    obj.D=-Am*Dpx-Ap*Dmx-Bm*Dpy-Bp*Dmy-Cm*Dpz-Cp*Dmz-obj.Eevaluated;
                otherwise
                    D1_x = kr(I_n, ops_x.D1, I_y,I_z);
                    D1_y = kr(I_n, I_x, ops_y.D1,I_z);
                    D1_z = kr(I_n, I_x, I_y,ops_z.D1);
                    obj.D=-obj.Aevaluated*D1_x-obj.Bevaluated*D1_y-obj.Cevaluated*D1_z-obj.Eevaluated;
            end
        end
        
        % Closure functions return the opertors applied to the own doamin to close the boundary
        % Penalty functions return the opertors to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of boundary condition if there are several.
        %       data                is a function returning the data that should be applied at the boundary.
        function [closure, penalty] = boundary_condition(obj,boundary,type,L)
            default_arg('type','char');
            BM=boundary_matrices(obj,boundary);
            
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
        
        function [ret] = evaluateCoefficientMatrix(obj, mat, X, Y, Z)
            params = obj.params;
            side = max(length(X),length(Y));
            if isa(mat,'function_handle')
                [rows,cols] = size(mat(params,0,0,0));
                matVec = mat(params,X',Y',Z');
                matVec=sparse(matVec);
            else
                matVec = mat;
                [rows,cols]=size(matVec);
                side=max(length(X),length(Y));
                cols=cols/side;
            end
            ret=cell(rows,cols);
            
            for ii=1:rows
                for jj=1:cols
                    ret{ii,jj}=diag(matVec(ii,(jj-1)*side+1:jj*side));
                end
            end
            ret=cell2mat(ret);
        end
        
        
        function [BM]=boundary_matrices(obj,boundary)
            params=obj.params;
            
            switch boundary
                case {'w','W','west'}
                    BM.e_=obj.e_w;
                    mat=obj.A;
                    BM.boundpos='l';
                    BM.Hi=obj.Hxi;
                    [BM.V,BM.Vi,BM.D,signVec]=obj.matrixDiag(mat,obj.X(1),obj.Yx,obj.Zx);
                    BM.side=length(obj.Yx);
                case {'e','E','east'}
                    BM.e_=obj.e_e;
                    mat=obj.A;
                    BM.boundpos='r';
                    BM.Hi=obj.Hxi;
                    [BM.V,BM.Vi,BM.D,signVec]=obj.matrixDiag(mat,obj.X(end),obj.Yx,obj.Zx);
                    BM.side=length(obj.Yx);
                case {'s','S','south'}
                    BM.e_=obj.e_s;
                    mat=obj.B;
                    BM.boundpos='l';
                    BM.Hi=obj.Hyi;
                    [BM.V,BM.Vi,BM.D,signVec]=obj.matrixDiag(mat,obj.Xy,obj.Y(1),obj.Zy);
                    BM.side=length(obj.Xy);
                case {'n','N','north'}
                    BM.e_=obj.e_n;
                    mat=obj.B;
                    BM.boundpos='r';
                    BM.Hi=obj.Hyi;
                    [BM.V,BM.Vi,BM.D,signVec]=obj.matrixDiag(mat,obj.Xy,obj.Y(end),obj.Zy);
                    BM.side=length(obj.Xy);
                case{'b','B','Bottom'}
                    BM.e_=obj.e_b;
                    mat=obj.C;
                    BM.boundpos='l';
                    BM.Hi=obj.Hzi;
                    [BM.V,BM.Vi,BM.D,signVec]=obj.matrixDiag(mat,obj.Xz,obj.Yz,obj.Z(1));
                    BM.side=length(obj.Xz);
                case{'t','T','Top'}
                    BM.e_=obj.e_t;
                    mat=obj.C;
                    BM.boundpos='r';
                    BM.Hi=obj.Hzi;
                    [BM.V,BM.Vi,BM.D,signVec]=obj.matrixDiag(mat,obj.Xz,obj.Yz,obj.Z(end));
                    BM.side=length(obj.Xz);
            end
            
            BM.pos=signVec(1); BM.zeroval=signVec(2); BM.neg=signVec(3);
        end
        
        
        function [closure, penalty]=boundary_condition_char(obj,BM)
            side = BM.side;
            pos = BM.pos;
            neg = BM.neg;
            zeroval=BM.zeroval;
            V = BM.V;
            Vi = BM.Vi;
            Hi=BM.Hi;
            D=BM.D;
            e_=BM.e_;
            
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
            zeroval=BM.zeroval;
            V = BM.V;
            Vi = BM.Vi;
            Hi=BM.Hi;
            D=BM.D;
            e_=BM.e_;
            switch boundary
                case {'w','W','west'}
                    L = obj.evaluateCoefficientMatrix(L,obj.x(1),obj.Yx,obj.Zx);
                case {'e','E','east'}
                    L = obj.evaluateCoefficientMatrix(L,obj.x(end),obj.Yx,obj.Zx);
                case {'s','S','south'}
                    L = obj.evaluateCoefficientMatrix(L,obj.Xy,obj.y(1),obj.Zy);
                case {'n','N','north'}
                    L = obj.evaluateCoefficientMatrix(L,obj.Xy,obj.y(end),obj.Zy);
                case {'b','B','bottom'}
                    L = obj.evaluateCoefficientMatrix(L,obj.Xz,obj.Yz,obj.z(1));
                case {'t','T','top'}
                    L = obj.evaluateCoefficientMatrix(L,obj.Xz,obj.Yz,obj.z(end));
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
        
        
        function [V,Vi, D,signVec]=matrixDiag(obj,mat,x,y,z)
            params = obj.params;
            syms xs ys zs
            [V, D] = eig(mat(params,xs,ys,zs));
            Vi=inv(V);
            xs = x;
            ys = y;
            zs = z;
            
            
            side = max(length(x),length(y));
            Dret = zeros(obj.n,side*obj.n);
            Vret = zeros(obj.n,side*obj.n);
            Viret= zeros(obj.n,side*obj.n);
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
            V = obj.evaluateCoefficientMatrix(V,x,y,z);
            Vi= obj.evaluateCoefficientMatrix(Vi,x,y,z);
            D = obj.evaluateCoefficientMatrix(D,x,y,z);
            DD = diag(D);
            
            poseig = (DD>0);
            zeroeig = (DD==0);
            negeig = (DD<0);
            
            D = diag([DD(poseig); DD(zeroeig); DD(negeig)]);
            V = [V(:,poseig) V(:,zeroeig) V(:,negeig)];
            Vi= [Vi(poseig,:); Vi(zeroeig,:); Vi(negeig,:)];
            signVec = [sum(poseig),sum(zeroeig),sum(negeig)];
        end
    end
end