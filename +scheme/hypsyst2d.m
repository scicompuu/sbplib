classdef hypsyst2d < scheme.Scheme
    properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing
        x,y % Grid
        X,Y % Values of x and y for each grid point
        order % Order accuracy for the approximation
        
        D % non-stabalized scheme operator
        A, B, E
        
        H % Discrete norm
        Hi
        H_x, H_y % Norms in the x and y directions
        Hx,Hy % Kroneckerd norms. 1'*Hx*v corresponds to integration in the x dir.
        I_x,I_y
        e_w, e_e, e_s, e_n
        params %parameters for the coeficient matrices
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
            
            ops_x = sbp.D2Standard(m_x,xlim,order);
            ops_y = sbp.D2Standard(m_y,ylim,order);
            
            obj.x=ops_x.x;
            obj.y=ops_y.y;
            
            obj.X = kr(x,ones(m_y,1));
            obj.Y = kr(ones(m_x,1),y);
            
            I_x = speye(m_x);
            I_y = speye(m_y);
            I_n=  eye(4);
            
            
            D1_x = kr(kr(I_n,ops_x.D1),I_y);
            obj.Hi_x= kr(kr(I_n,ops_x.HI),I_y);
            D1_y=kr(I_n,kr(I_x,ops_y.D1));
            obj.Hi_y=kr(I_n,kr(I_x,ops_y.HI));
            
            obj.e_w=kr(I_n,kr(ops_x.e_l,I_y));
            obj.e_e=kr(I_n,kr(ops_x.e_r,I_y));
            obj.e_s=kr(I_n,kr(I_x,ops_y.e_l));
            obj.e_n=kr(I_n,kr(I_x,ops_y.e_r));
            
            obj.m=m;
            obj.h=[ops_x.h ops_y.h];
            obj.order=order;
            obj.params=params;
            
            obj.A=matrixBuild(obj,matrices.A);
            obj.B=matrixBuild(obj,matrices.B);
            obj.E=matrixBuild(obj,matrices.E);
            
            obj.D=-obj.A*D1_x-obj.B*D1_y-E;
            
        end
        % Closure functions return the opertors applied to the own doamin to close the boundary
        % Penalty functions return the opertors to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       type                is a string specifying the type of boundary condition if there are several.
        %       data                is a function returning the data that should be applied at the boundary.
        %       neighbour_scheme    is an instance of Scheme that should be interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to interface to.
        function [closure, penalty] = boundary_condition(obj,boundary,type,data)
            default_arg('type','neumann');
            default_arg('data',0);
            
            switch type
                case{c,'char'}
                    [tau,e_,Hi,CHM]=GetBoundarydata(obj,boundary,type);
                    closure =Hi*e_*tau*CHM*e_';
                    penalty =Hi*e_*tau*CHM*data;
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
        
    end
    
    methods(Static)
        % Calculates the matrcis need for the inteface coupling between boundary bound_u of scheme schm_u
        % and bound_v of scheme schm_v.
        %   [uu, uv, vv, vu] = inteface_couplong(A,'r',B,'l')
        function [uu, uv, vv, vu] = interface_coupling(schm_u,bound_u,schm_v,bound_v)
            [uu,uv] = schm_u.interface(bound_u,schm_v,bound_v);
            [vv,vu] = schm_v.interface(bound_v,schm_u,bound_u);
        end
        
        function [ret]=matrixBuild(obj,mat)
            %extra info for coordinate transfomration mult my y_ny  and
            %x,ny osv...
            params=obj.params;
            X=obj.X;
            Y=obj.Y;
            
            if isa(mat,'function_handle')
                matVec=mat(params,x,y);
                side=length(x);
            else
                matVec=mat;
                side=max(size(mat))/4;
            end
            
            
            ret=[diag(matVec(1,(1-1)*side+1:1*side)) diag(matVec(1,(2-1)*side+1:2*side)) diag(matVec(1,(3-1)*side+1:3*side)) diag(matVec(1,(4-1)*side+1:4*side))
                diag(matVec(2,(1-1)*side+1:1*side)) diag(matVec(2,(2-1)*side+1:2*side)) diag(matVec(2,(3-1)*side+1:3*side)) diag(matVec(2,(4-1)*side+1:4*side))
                diag(matVec(3,(1-1)*side+1:1*side)) diag(matVec(3,(2-1)*side+1:2*side)) diag(matVec(3,(3-1)*side+1:3*side)) diag(matVec(3,(4-1)*side+1:4*side))
                diag(matVec(4,(1-1)*side+1:1*side)) diag(matVec(4,(2-1)*side+1:2*side)) diag(matVec(4,(3-1)*side+1:3*side)) diag(matVec(4,(4-1)*side+1:4*side))];
        end
        
        function [tau,e_,Hi, CHM]=GetBoundarydata(obj,boundary)
            params=obj.params;
            x=obj.x;
            y=obj.y;
            
            side=max(length(x),length(y));
            
            
            switch boundary
                case {'w','W','west'}
                    e_=obj.e_w;
                    mat=obj.A;
                    [V,D]=matrixDiag(mat,params,x(1),y);
                    Hi=obj.Hx;
                    tau=zeros(4*side,pos*side);
                    tau(1:pos*side,:)=-abs(D(1:pos*side,1:pos*side));
                    CHM=V*((D+abs(D))/2)*V';
                case {'e','E','east'}
                    e_=obj.e_e;
                    mat=obj.A;
                    [V,D]=matrixDiag(mat,params,x(end),y);
                    Hi=obj.Hy;
                    tau=zeros(4*side,(4-pos));
                    tau(pos*side+1:4*side)=-abs(D(pos*side+1:4*side,pos*side+1:4*side));
                    CHM=V*((D-abs(D))/2)*V';
                case {'s','S','south'}
                    e_=obj.e_s;
                    mat=obj.B;
                    [V,D]=matrixDiag(mat,params,x,y(1));
                    Hi=obj.Hx;
                    tau=zeros(4*side,pos*side);
                    tau(1:pos*side,:)=-abs(D(1:pos*side,1:pos*side));
                    CHM=V*((D+abs(D))/2)*V';
                case {'n','N','north'}
                    e_=obj.e_n;
                    mat=obj.B;
                    [V,D]=matrixDiag(mat,params,x,y(end));
                    Hi=obj.Hy;
                    tau=zeros(4*side,(4-pos));
                    tau(pos*side+1:4*side)=-abs(D(pos*side+1:4*side,pos*side+1:4*side));
                    CHM=V*((D-abs(D))/2)*V';
            end
            
            tau=V*tau*V';
            
        end
        
        function [V, D,pos]=matrixDiag(mat,params,x,y)
            syms xs ys;
            [V, D]=eig(mat(params,xs,ys));
            xs=1;ys=1;
            DD=eval(diag(D));
            
            pos=find(DD>=0); %Now zero eigenvalues are calculated as possitive, Maybe it should not????
            neg=find(DD<0);
            syms xs ys
            DD=diag(D);
            
            D=diag([DD(pos); DD(neg)]);
            V=[V(:,pos) V(:,neg)];
            
            xs=x; ys=y;
            
            side=max(length(x),length(y));
            Dret=zeros(4,side*4);
            Vret=zeros(4,side*4);
            for ii=1:4
                for jj=1:4
                    Dret(jj,(ii-1)*side+1:side*ii)=eval(D(jj,ii));
                    Vret(jj,(ii-1)*side+1:side*ii)=eval(V(jj,ii));
                end
            end
            V=sparse(normc(Vret));
            D=sparse(Dret);
            
            
            V=matrixBuild([],[],[],V);
            D=matrixBuild([],[],[],D);
            pos=legth(pos);
        end
    end
end