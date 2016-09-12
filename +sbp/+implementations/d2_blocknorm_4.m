function [H, HI, D1, D2, e_1, e_m, M, Q, S_1, S_m] = d2_blocknorm_4(m,h)
    
    BP = 4;
    if(m<2*BP)
        error(['Operator requires at least ' num2str(2*BP) ' grid points']);
    end

    H_U=[0.751e3 / 0.3456e4 0.661e3 / 0.3456e4 -0.515e3 / 0.3456e4 0.5e1 / 0.128e3; 0.661e3 / 0.3456e4 0.1405e4 / 0.1152e4 -0.3e1 / 0.128e3 0.29e2 / 0.3456e4; -0.515e3 / 0.3456e4 -0.3e1 / 0.128e3 0.989e3 / 0.1152e4 0.149e3 / 0.3456e4; 0.5e1 / 0.128e3 0.29e2 / 0.3456e4 0.149e3 / 0.3456e4 0.3407e4 / 0.3456e4;];


    H=speye(m);
    H(1:4,1:4)=H_U;
    H(m-3:m,m-3:m)=rot90( H_U(1:4,1:4) ,2 );
    H=H*h;
    HI=inv(H);

    e=ones(m,1);
    Q=spdiags([e -8*e 0*e 8*e -e], -2:2, m, m)/12;
%     Q=(-1/12*diag(ones(m-2,1),2)+8/12*diag(ones(m-1,1),1)-8/12*diag(ones(m-1,1),-1)+1/12*diag(ones(m-2,1),-2));

    Q_U = [-0.1e1 / 0.2e1 0.55e2 / 0.72e2 -0.47e2 / 0.144e3 0.1e1 / 0.16e2; -0.55e2 / 0.72e2 0 0.43e2 / 0.48e2 -0.19e2 / 0.144e3; 0.47e2 / 0.144e3 -0.43e2 / 0.48e2 0 0.47e2 / 0.72e2; -0.1e1 / 0.16e2 0.19e2 / 0.144e3 -0.47e2 / 0.72e2 0;];

    Q(1:4,1:4)=Q_U;
    Q(m-3:m,m-3:m)=rot90( -Q_U(1:4,1:4) ,2 );

    D1=H\Q;

    M_U=[0.359e3 / 0.288e3 -0.443e3 / 0.288e3 0.97e2 / 0.288e3 -0.13e2 / 0.288e3; -0.51e2 / 0.32e2 0.325e3 / 0.96e2 -0.191e3 / 0.96e2 0.19e2 / 0.96e2; 0.43e2 / 0.96e2 -0.69e2 / 0.32e2 0.293e3 / 0.96e2 -0.137e3 / 0.96e2; -0.29e2 / 0.288e3 0.89e2 / 0.288e3 -0.427e3 / 0.288e3 0.727e3 / 0.288e3;];



%     M=-(-1/12*diag(ones(m-2,1),2)+16/12*diag(ones(m-1,1),1)+16/12*diag(ones(m-1,1),-1)-1/12*diag(ones(m-2,1),-2)-30/12*diag(ones(m,1),0));
    M=-spdiags([-e 16*e -30*e 16*e -e], -2:2, m, m)/12;

    M(1:4,1:4)=M_U;

    M(m-3:m,m-3:m)=rot90(  M_U ,2 );
    M=M/h;

    DS_U=[0.25e2 / 0.12e2 -4 3 -0.4e1 / 0.3e1 0.1e1 / 0.4e1;];
    DS=sparse(m,m);
    DS(1,1:5)=DS_U;
    DS(m,m-4:m)=fliplr(DS_U);
    DS=DS/h;

    D2=H\(-M+DS);

%     d3=[-1 3 -3 1];
%     t3=sum(abs(d3));
%     DD_3(1:1,1:4)=[d3];
%     DD_3(m:m,m-3:m)=[d3];

    % This works for wave eq.
    % For studs interface in 1D no AD is needed.
%     ADD=1*h/(t3)*DD_3'*DD_3;

    e_1 = sparse(m,1);
    e_1(1)= 1;
    e_m = sparse(m,1);
    e_m(end)= 1;
    S_1 = -DS(1,:)';
    S_m =  DS(end,:)';

    Q = H*D1-(-(e_1*e_1') + (e_m*e_m'));
    M = -(H*D2-(-e_1*S_1' + e_m*S_m'));
end