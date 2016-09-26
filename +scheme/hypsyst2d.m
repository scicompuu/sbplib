classdef hypsyst2d < scheme.Scheme
    properties
        m % Number of points in each direction, possibly a vector
        n %size of system
        h % Grid spacing
        x,y % Grid
        X,Y % Values of x and y for each grid point
        order % Order accuracy for the approximation
        
        D % non-stabalized scheme operator
        A, B, E
        
        H % Discrete norm
        % Norms in the x and y directions
        Hxi,Hyi % Kroneckerd norms. 1'*Hx*v corresponds to integration in the x dir.
        I_x,I_y, I_N
        e_w, e_e, e_s, e_n
        params %parameters for the coeficient matrices
        matrices
    end
    
    
    methods
        function obj = hypsyst2d(m,lim,order,matrices,params)
            
            xlim = lim{1};
            ylim = lim{2};
            
            if length(m) == 1
                m = [m m];
            end
            
            m_x = m(1);
            m_y = m(2);
            obj.params=params;
            
            obj.matrices=matrices;
            
            ops_x = sbp.D2Standard(m_x,xlim,order);
            ops_y = sbp.D2Standard(m_y,ylim,order);
            
            obj.x=ops_x.x;
            obj.y=ops_y.x;
            
            obj.X = kr(obj.x,ones(m_y,1));
            obj.Y = kr(ones(m_x,1),obj.y);
            
            obj.A=obj.matrixBuild(matrices.A);
            obj.B=obj.matrixBuild(matrices.B);
            obj.E=obj.matrixBuild(matrices.E);
            
            obj.n=length(matrices.A(obj.params,0,0));
            
            I_n=  eye(obj.n);
            I_x = speye(m_x); obj.I_x=I_x;
            I_y = speye(m_y); obj.I_y=I_y;
            
            
            D1_x = kr(kr(I_n,ops_x.D1),I_y);
            obj.Hxi= kr(kr(I_n,ops_x.HI),I_y);
            D1_y=kr(I_n,kr(I_x,ops_y.D1));
            obj.Hyi=kr(I_n,kr(I_x,ops_y.HI));
            
            obj.e_w=kr(I_n,kr(ops_x.e_l,I_y));
            obj.e_e=kr(I_n,kr(ops_x.e_r,I_y));
            obj.e_s=kr(I_n,kr(I_x,ops_y.e_l));
            obj.e_n=kr(I_n,kr(I_x,ops_y.e_r));
            
            obj.m=m;
            obj.h=[ops_x.h ops_y.h];
            obj.order=order;
            
            obj.D=-obj.A*D1_x-obj.B*D1_y-obj.E;
            
        end
        
        % Closure functions return the opertors applied to the own doamin to close the boundary
        % Penalty functions return the opertors to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of boundary condition if there are several.
        %       data                is a function returning the data that should be applied at the boundary.
        %       neighbour_scheme    is an instance of Scheme that should be interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to interface to.
        function [closure, penalty] = boundary_condition(obj,boundary,type,L)
            default_arg('type','char');
            switch type
                case{'c','char'}
                    [closure,penalty]=GetBoundarydata_char(obj,boundary);
                case{'general'}
                    [closure,penalty]=GeneralBoundaryCond(obj,boundary,L);
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
        
        function [ret]=matrixBuild(obj,mat,X,Y)
            params=obj.params;
            default_arg('X',obj.X);
            default_arg('Y',obj.Y)
            
            if isa(mat,'function_handle')
                [rows,cols]=size(mat(params,0,0));
                matVec=mat(params,X',Y');
                matVec=sparse(matVec);
                side=max(length(X),length(Y));
            else
                matVec=mat;
                [rows,cols]=size(matVec);
                side=max(length(X),length(Y));
                cols=cols/side;
            end
            ret=kron(ones(rows,cols),speye(side));
            
            for ii=1:rows
                for jj=1:cols
                    ret((ii-1)*side+1:ii*side,(jj-1)*side+1:jj*side)=diag(matVec(ii,(jj-1)*side+1:jj*side));
                end
            end
        end
        
        
        function [closure, penalty]=GetBoundarydata_char(obj,boundary)
            params=obj.params;
            x=obj.x; y=obj.y;
            side=max(length(x),length(y));
            
            switch boundary
                case {'w','W','west'}
                    e_=obj.e_w;
                    mat=obj.matrices.A;
                    boundPos='l';
                    Hi=obj.Hxi;
                    [V,Vi,D,signVec]=obj.matrixDiag(mat,x(1),y);
                case {'e','E','east'}
                    e_=obj.e_e;
                    mat=obj.matrices.A;
                    boundPos='r';
                    Hi=obj.Hxi;
                    [V,Vi,D,signVec]=obj.matrixDiag(mat,x(end),y);
                case {'s','S','south'}
                    e_=obj.e_s;
                    mat=obj.matrices.B;
                    boundPos='l';
                    Hi=obj.Hxi;
                    [V,Vi,D,signVec]=obj.matrixDiag(mat,x,y(1));
                case {'n','N','north'}
                    e_=obj.e_n;
                    mat=obj.matrices.B;
                    boundPos='r';
                    Hi=obj.Hxi;
                    [V,Vi,D,signVec]=obj.matrixDiag(mat,x,y(end));
            end
            
            pos=signVec(1); zeroval=signVec(2); neg=signVec(3);
            
            switch boundPos
                case {'l'}
                    tau=sparse(obj.n*side,pos*side);
                    Vi_plus=Vi(1:pos*side,:);
                    tau(1:pos*side,:)=-abs(D(1:pos*side,1:pos*side));
                    closure=Hi*e_*V*tau*Vi_plus*e_';
                    penalty=-Hi*e_*V*tau*Vi_plus;
                case {'r'}
                    tau=sparse(obj.n*side,neg*side);
                    tau((pos+zeroval)*side+1:obj.n*side,:)=-abs(D((pos+zeroval)*side+1:obj.n*side,(pos+zeroval)*side+1:obj.n*side));
                    Vi_minus=Vi((pos+zeroval)*side+1:obj.n*side,:);
                    closure=Hi*e_*V*tau*Vi_minus*e_';
                    penalty=-Hi*e_*V*tau*Vi_minus;
            end
        end
        
        
        function [closure,penalty]=GeneralBoundaryCond(obj,boundary,L)
            params=obj.params;
            x=obj.x; y=obj.y;
            side=max(length(x),length(y));
            
            switch boundary
                case {'w','W','west'}
                    e_=obj.e_w;
                    mat=obj.matrices.A;
                    boundPos='l';
                    Hi=obj.Hxi;
                    [V,Vi,D,signVec]=obj.matrixDiag(mat,x(1),y);
                    L=obj.matrixBuild(L,x(1),y);
                case {'e','E','east'}
                    e_=obj.e_e;
                    mat=obj.matrices.A;
                    boundPos='r';
                    Hi=obj.Hxi;
                    [V,Vi,D,signVec]=obj.matrixDiag(mat,x(end),y);
                    L=obj.matrixBuild(L,x(end),y);
                case {'s','S','south'}
                    e_=obj.e_s;
                    mat=obj.matrices.B;
                    boundPos='l';
                    Hi=obj.Hxi;
                    [V,Vi,D,signVec]=obj.matrixDiag(mat,x,y(1));
                    L=obj.matrixBuild(L,x,y(1));
                case {'n','N','north'}
                    e_=obj.e_n;
                    mat=obj.matrices.B;
                    boundPos='r';
                    Hi=obj.Hxi;
                    [V,Vi,D,signVec]=obj.matrixDiag(mat,x,y(end));
                    L=obj.matrixBuild(L,x,y(end));
            end
            
            pos=signVec(1); zeroval=signVec(2); neg=signVec(3);
            
            switch boundPos
                case {'l'}
                    tau=sparse(obj.n*side,pos*side);
                    Vi_plus=Vi(1:pos*side,:);
                    Vi_minus=Vi(pos*side+1:obj.n*side,:);
                    V_plus=V(:,1:pos*side);
                    V_minus=V(:,(pos+zeroval)*side+1:obj.n*side);
                    
                    tau(1:pos*side,:)=-abs(D(1:pos*side,1:pos*side));
                    R=-inv(L*V_plus)*(L*V_minus);
                    closure=Hi*e_*V*tau*(Vi_plus-R*Vi_minus)*e_';
                    penalty=-Hi*e_*V*tau*inv(L*V_plus)*L;
                case {'r'}
                    tau=sparse(obj.n*side,neg*side);
                    tau((pos+zeroval)*side+1:obj.n*side,:)=-abs(D((pos+zeroval)*side+1:obj.n*side,(pos+zeroval)*side+1:obj.n*side));
                    Vi_plus=Vi(1:pos*side,:);
                    Vi_minus=Vi((pos+zeroval)*side+1:obj.n*side,:);
                    
                    V_plus=V(:,1:pos*side);
                    V_minus=V(:,(pos+zeroval)*side+1:obj.n*side);
                    R=-inv(L*V_minus)*(L*V_plus);
                    closure=Hi*e_*V*tau*(Vi_minus-R*Vi_plus)*e_';
                    penalty=-Hi*e_*V*tau*inv(L*V_minus)*L;
            end
        end
        
        
        function [V,Vi, D,signVec]=matrixDiag(obj,mat,x,y)
            params=obj.params;
            syms xs ys;
            [V, D]=eig(mat(params,xs,ys));
            xs=1;ys=1;
            DD=eval(diag(D));
            
            poseig=find(DD>0);
            zeroeig=find(DD==0);
            negeig=find(DD<0);
            syms xs ys
            DD=diag(D);
            
            D=diag([DD(poseig);DD(zeroeig); DD(negeig)]);
            V=[V(:,poseig) V(:,zeroeig) V(:,negeig)];
            xs=x; ys=y;
            
            side=max(length(x),length(y));
            Dret=zeros(obj.n,side*obj.n);
            Vret=zeros(obj.n,side*obj.n);
            for ii=1:obj.n
                for jj=1:obj.n
                    Dret(jj,(ii-1)*side+1:side*ii)=eval(D(jj,ii));
                    Vret(jj,(ii-1)*side+1:side*ii)=eval(V(jj,ii));
                end
            end
            
            D=sparse(Dret);
            V=sparse(normc(Vret));
            V=obj.matrixBuild(V,x,y);
            D=obj.matrixBuild(D,x,y);
            Vi=inv(V);
            signVec=[length(poseig),length(zeroeig),length(negeig)];
        end
        
    end
    
    methods(Static)
        % Calculates the matrcis need for the inteface coupling between boundary bound_u of scheme schm_u
        % and bound_v of scheme schm_v.
        %   [uu, uv, vv, vu] = inteface_couplong(A,'r',B,'l')
        function [uu, uv, vv, vu] = interface_coupling(schm_u,bound_u,schm_v,bound_v)
            [uu,uv] = schm_u.interface(bound_u,schm_v,bound_v);
            [vv,vu] = schm_v.interface(bound_v,schm_u,bound_u);
        end
        
        
    end
end