function [H, HI, D1, D2, e_l, e_r, d1_l, d1_r] = d2_variable_2(m,h)

    BP = 1;
    if(m<2*BP)
        error(['Operator requires at least ' num2str(2*BP) ' grid points']);
    end

    % Norm
    Hv = ones(m,1);
    Hv(1) = 1/2;
    Hv(m:m) = 1/2;
    Hv = h*Hv;
    H = spdiag(Hv, 0);
    HI = spdiag(1./Hv, 0);


    % Boundary operators
    e_l = sparse(m,1);
    e_l(1) = 1;
    e_r = rot90(e_l, 2);

    d1_l = sparse(m,1);
    d1_l(1:3) = 1/h*[-3/2 2 -1/2];
    d1_r = -rot90(d1_l, 2);

    % D1 operator
    diags   = -1:1;
    stencil = [-1/2 0 1/2];
    D1 = stripeMatrix(stencil, diags, m);
    
    D1(1,1)=-1;D1(1,2)=1;D1(m,m-1)=-1;D1(m,m)=1;
    D1(m,m-1)=-1;D1(m,m)=1;
    D1=D1/h;
    %Q=H*D1 + 1/2*(e_1*e_1') - 1/2*(e_m*e_m');


    M=sparse(m,m);

    scheme_width = 3;
    scheme_radius = (scheme_width-1)/2;
    r = (1+scheme_radius):(m-scheme_radius);

    function D2 = D2_fun(c)

        Mm1 = -c(r-1)/2 - c(r)/2;
        M0  =  c(r-1)/2 + c(r)   + c(r+1)/2;
        Mp1 =            -c(r)/2 - c(r+1)/2;

        M(r,:) = spdiags([Mm1 M0 Mp1],0:2*scheme_radius,length(r),m);


        M(1:2,1:2)=[c(1)/2 + c(2)/2 -c(1)/2 - c(2)/2; -c(1)/2 - c(2)/2 c(1)/2 + c(2) + c(3)/2;];
        M(m-1:m,m-1:m)=[c(m-2)/2 + c(m-1) + c(m)/2 -c(m-1)/2 - c(m)/2; -c(m-1)/2 - c(m)/2 c(m-1)/2 + c(m)/2;];
        M=M/h;

        D2=HI*(-M-c(1)*e_l*d1_l'+c(m)*e_r*d1_r');
    end
    D2 = @D2_fun;
end