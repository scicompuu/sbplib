classdef Hypsyst2dCurve < scheme.Scheme
    properties
        m % Number of points in each direction, possibly a vector
        n %size of system
        h % Grid spacing
        X,Y % Values of x and y for each grid point
        
        J, Ji
        xi,eta
        Xi,Eta
        
        A,B
        X_eta, Y_eta 
        X_xi,Y_xi
        order % Order accuracy for the approximation

        D % non-stabalized scheme operator
        Ahat, Bhat, E
    
        H % Discrete norm
        % Norms in the x and y directions
        Hxii,Hetai % Kroneckerd norms. 1'*Hx*v corresponds to integration in the x dir.
        I_xi,I_eta, I_N, onesN
        e_w, e_e, e_s, e_n
        index_w, index_e,index_s,index_n
        params %parameters for the coeficient matrice
    end


    methods
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
            
            obj.index_w=1:m_xi;
            obj.index_e=(m_tot-m_xi+1):m_tot;
            obj.index_s=1:m_xi:(m_tot-m_xi+1);
            obj.index_n=(m_xi):m_xi:m_tot;

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
            obj.e_n = kr(I_n, I_xi, ops_eta.e_r);
            
            [X,Y] = ti.map(obj.xi,obj.eta);
                   
            [x_xi,x_eta] = gridDerivatives(X,ops_xi.D1,ops_eta.D1);
            [y_xi,y_eta] = gridDerivatives(Y,ops_xi.D1, ops_eta.D1);
                    
            obj.X=reshape(X,m_xi*m_eta,1);
            obj.Y=reshape(Y,m_xi*m_eta,1);
            obj.X_xi=reshape(x_xi,m_xi*m_eta,1);
            obj.Y_xi=reshape(y_xi,m_xi*m_eta,1);
            obj.X_eta=reshape(x_eta,m_xi*m_eta,1);
            obj.Y_eta=reshape(y_eta,m_xi*m_eta,1);
           
            Ahat_evaluated = obj.evaluateCoefficientMatrix(obj.Ahat, obj.X, obj.Y,obj.X_eta,obj.Y_eta);
            Bhat_evaluated = obj.evaluateCoefficientMatrix(obj.Bhat, obj.X, obj.Y,obj.X_xi,obj.Y_xi);
            E_evaluated = obj.evaluateCoefficientMatrix(obj.E, obj.X, obj.Y,[],[]);

            obj.m=m;
            obj.h=[ops_xi.h ops_eta.h];
            obj.order=order;
            obj.J=x_xi.*y_eta - x_eta.*y_xi;  
            obj.Ji =kr(I_n,spdiags(1./obj.J(:),0,m_xi*m_eta,m_xi*m_eta));

            obj.D=obj.Ji*(-Ahat_evaluated*D1_xi-Bhat_evaluated*D1_eta)-E_evaluated;

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
                    [closure,penalty]=boundary_condition_char(obj,boundary);
                case{'general'}
                    [closure,penalty]=boundary_condition_general(obj,boundary,L);
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

        function [ret] = evaluateCoefficientMatrix(obj, mat, X, Y,x_,y_)
            params=obj.params;

            if isa(mat,'function_handle')
                [rows,cols]=size(mat(params,0,0,0,0));
                x_=kr(x_,obj.onesN);
                y_=kr(y_,obj.onesN);
                matVec=mat(params,X',Y',x_',y_');
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


        function [closure, penalty]=boundary_condition_char(obj,boundary)
            params=obj.params;
            xi=obj.xi; eta=obj.eta;
            side=max(length(xi),length(eta));

            switch boundary
                case {'w','W','west'}
                    e_=obj.e_w;
                    mat=obj.Ahat;
                    boundPos='l';
                    Hi=obj.Hxii;
                    [V,Vi,D,signVec]=obj.matrixDiag(mat,xi(1),eta,obj.X_eta(obj.index_w),obj.Y_eta(obj.index_w));
                case {'e','E','east'}
                    e_=obj.e_e;
                    mat=obj.Ahat;
                    boundPos='r';
                    Hi=obj.Hxii;
                    [V,Vi,D,signVec]=obj.matrixDiag(mat,xi(end),eta,obj.X_eta(obj.index_e),obj.Y_eta(obj.index_e));
                case {'s','S','south'}
                    e_=obj.e_s;
                    mat=obj.Bhat;
                    boundPos='l';
                    Hi=obj.Hetai;
                    [V,Vi,D,signVec]=obj.matrixDiag(mat,xi,eta(1),obj.X_xi(obj.index_s),obj.Y_xi(obj.index_s));
                case {'n','N','north'}
                    e_=obj.e_n;
                    mat=obj.Bhat;
                    boundPos='r';
                    Hi=obj.Hetai;
                    [V,Vi,D,signVec]=obj.matrixDiag(mat,xi,eta(end),obj.X_xi(obj.index_n),obj.Y_xi(obj.index_n));
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


        function [closure,penalty]=boundary_condition_general(obj,boundary,L)
            params=obj.params;
            xi=obj.xi; eta=obj.eta;
            side=max(length(xi),length(eta));

            switch boundary
                case {'w','W','west'}
                    e_=obj.e_w;
                    mat=obj.Ahat;
                    boundPos='l';
                    Hi=obj.Hxii;
                    [V,Vi,D,signVec]=obj.matrixDiag(mat,xi(1),eta,obj.x_eta,obj.y_eta);
                    L=obj.evaluateCoefficientMatrix(L,xi(1),eta);
                case {'e','E','east'}
                    e_=obj.e_e;
                    mat=obj.Ahat; 
                    boundPos='r';
                    Hi=obj.Hxii;
                    [V,Vi,D,signVec]=obj.matrixDiag(mat,xi(end),eta,obj.x_eta,obj.y_eta);
                    L=obj.evaluateCoefficientMatrix(L,xi(end),eta,[],[]);
                case {'s','S','south'}
                    e_=obj.e_s;
                    mat=obj.Bhat;
                    boundPos='l';
                    Hi=obj.Hetai;
                    [V,Vi,D,signVec]=obj.matrixDiag(mat,xi,eta(1),obj.x_xi,obj.y_xi);
                    L=obj.evaluateCoefficientMatrix(L,xi,eta(1));
                case {'n','N','north'}
                    e_=obj.e_n;
                    mat=obj.Bhat;
                    boundPos='r';
                    Hi=obj.Hetai;
                    [V,Vi,D,signVec]=obj.matrixDiag(mat,xi,eta(end),obj.x_xi,obj.y_xi);
                    L=obj.evaluateCoefficientMatrix(L,xi,eta(end));
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


        function [V,Vi, D,signVec]=matrixDiag(obj,mat,x,y,x_,y_)
            params=obj.params;
            syms xs ys 
            if(sum(abs(x_))~=0)
                syms xs_
            else
                xs_=0;
            end
            
            if(sum(abs(y_))~=0)
            syms ys_;
            else
                ys_=0;
            end
            
            [V, D]=eig(mat(params,xs,ys,xs_,ys_));
            xs=1;ys=1; xs_=x_(1); ys_=y_(1);
            DD=eval(diag(D));

            poseig=find(DD>0);
            zeroeig=find(DD==0);
            negeig=find(DD<0);
            syms xs ys xs_ ys_
            DD=diag(D);

            D=diag([DD(poseig);DD(zeroeig); DD(negeig)]);
            V=[V(:,poseig) V(:,zeroeig) V(:,negeig)];
            Vi=inv(V);
            xs=x; 
            ys=y;
            xs_=x_';
            ys_=y_';

            side=max(length(x),length(y));
            Dret=zeros(obj.n,side*obj.n);
            Vret=zeros(obj.n,side*obj.n);
            Viret=zeros(obj.n,side*obj.n);
            for ii=1:obj.n
                for jj=1:obj.n
                    Dret(jj,(ii-1)*side+1:side*ii)=eval(D(jj,ii));
                    Vret(jj,(ii-1)*side+1:side*ii)=eval(V(jj,ii));
                    Viret(jj,(ii-1)*side+1:side*ii)=eval(Vi(jj,ii));
                end
            end

            D=sparse(Dret);
            V=sparse(Vret);
            Vi=sparse(Viret);
            V=obj.evaluateCoefficientMatrix(V,x,y,x_,y_);
            D=obj.evaluateCoefficientMatrix(D,x,y,x_,y_);
            Vi=obj.evaluateCoefficientMatrix(Vi,x,y,x_,y_);
            signVec=[length(poseig),length(zeroeig),length(negeig)];
        end

    end
end