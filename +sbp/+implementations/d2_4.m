function [H, HI, D1, D2, e_1, e_m, M, Q, S_1, S_m] = d2_4(m,h)

    BP = 4;
    if(m<2*BP)
        error(['Operator requires at least ' num2str(2*BP) ' grid points']);
    end

    e_1=sparse(m,1);e_1(1)=1;
    e_m=sparse(m,1);e_m(m)=1;

    e=ones(m,1);
    H=spdiags(e,0,m,m);
    %H=speye(m,m);

    H(1:4,1:4)=diag([17/48 59/48 43/48 49/48]);
    H(m-3:m,m-3:m)=rot90(diag([17/48 59/48 43/48 49/48]),2);
    H=H*h;
    HI=inv(H);


    Q=spdiags([e -8*e 0*e 8*e -e], -2:2, m, m)/12;
    %Q=(-1/12*diag(ones(m-2,1),2)+8/12*diag(ones(m-1,1),1)-8/12*diag(ones(m-1,1),-1)+1/12*diag(ones(m-2,1),-2));
    Q_U = [
        0              0.59e2/0.96e2  -0.1e1/0.12e2  -0.1e1/0.32e2;
        -0.59e2/0.96e2 0              0.59e2/0.96e2  0;
        0.1e1/0.12e2   -0.59e2/0.96e2 0              0.59e2/0.96e2;
        0.1e1/0.32e2   0              -0.59e2/0.96e2 0;
    ];
    Q(1:4,1:4)=Q_U;
    Q(m-3:m,m-3:m)=rot90( -Q_U(1:4,1:4) ,2 );

    D1=H\(Q-1/2*(e_1*e_1')+1/2*(e_m*e_m')) ;

    M=-spdiags([-e 16*e -30*e 16*e -e], -2:2, m, m)/12;

    %M=-(-1/12*diag(ones(m-2,1),2)+16/12*diag(ones(m-1,1),1)+16/12*diag(ones(m-1,1),-1)-1/12*diag(ones(m-2,1),-2)-30/12*diag(ones(m,1),0));

    M_U=[
        0.9e1/0.8e1 -0.59e2/0.48e2 0.1e1/0.12e2 0.1e1/0.48e2;
        -0.59e2/0.48e2 0.59e2/0.24e2 -0.59e2/0.48e2 0;
        0.1e1/0.12e2 -0.59e2/0.48e2 0.55e2/0.24e2 -0.59e2/0.48e2;
        0.1e1/0.48e2 0 -0.59e2/0.48e2 0.59e2/0.24e2;
    ];
    M(1:4,1:4)=M_U;

    M(m-3:m,m-3:m)=rot90(  M_U ,2 );
    M=M/h;

    S_U=[-0.11e2/0.6e1 3 -0.3e1/0.2e1 0.1e1/0.3e1;]/h;
    S_1=sparse(sparse(1,m));
    S_1(1:4)=S_U;
    S_m=sparse(sparse(1,m));
    S_m(m-3:m)=fliplr(-S_U);


    D2=H\(-M-e_1*S_1+e_m*S_m);
    S_1 = S_1';
    S_m = S_m';
end