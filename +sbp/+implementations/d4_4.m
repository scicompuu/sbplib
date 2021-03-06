function [H, HI, D1, D2, D3, D4, e_1, e_m, M, ...
    M4,Q, Q3, S2_1, S2_m, S3_1, S3_m, S_1, S_m] = d4_4(m,h)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% 4:de ordn. SBP Finita differens         %%%
    %%% operatorer framtagna av Ken Mattsson    %%%
    %%%                                         %%%
    %%% 6 randpunkter, diagonal norm            %%%
    %%%                                         %%%
    %%% Datum: 2013-11-11                       %%%
    %%%                                         %%%
    %%%                                         %%%
    %%% H           (Normen)                    %%%
    %%% D1          (approx f?rsta derivatan)   %%%
    %%% D2          (approx andra derivatan)    %%%
    %%% D3          (approx tredje derivatan)   %%%
    %%% D2          (approx fj?rde derivatan)   %%%
    %%%                                         %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % M?ste ange antal punkter (m) och stegl?ngd (h)
    % Notera att dessa opetratorer ?r framtagna f?r att anv?ndas n?r
    % vi har 3de och 4de derivator i v?r PDE
    % I annat fall anv?nd de "traditionella" som har noggrannare
    % randsplutningar f?r D1 och D2

    % Vi b?rjar med normen. Notera att alla SBP operatorer delar samma norm,
    % vilket ?r n?dv?ndigt f?r stabilitet
    
    BP = 6;
    if(m<2*BP)
        error(['Operator requires at least ' num2str(2*BP) ' grid points']);
    end

    H=speye(m,m);
    H_U=[0.35809e5 / 0.100800e6 0 0 0 0 0; 0 0.13297e5 / 0.11200e5 0 0 0 0; 0 0 0.5701e4 / 0.5600e4 0 0 0; 0 0 0 0.45109e5 / 0.50400e5 0 0; 0 0 0 0 0.35191e5 / 0.33600e5 0; 0 0 0 0 0 0.33503e5 / 0.33600e5;];

    H(1:6,1:6)=H_U;
    H(m-5:m,m-5:m)=rot90(H_U,2);
    H=H*h;
    HI=inv(H);


    % First derivative SBP operator, 1st order accurate at first 6 boundary points

%     q2=-1/12;q1=8/12;
%     Q=q2*(diag(ones(m-2,1),2) - diag(ones(m-2,1),-2))+q1*(diag(ones(m-1,1),1)-diag(ones(m-1,1),-1));
    e=ones(m,1);
    Q=spdiags([e -8*e 0*e 8*e -e], -2:2, m, m)/12;

    %Q=(-1/12*diag(ones(m-2,1),2)+8/12*diag(ones(m-1,1),1)-8/12*diag(ones(m-1,1),-1)+1/12*diag(ones(m-2,1),-2));

    Q_U = [0 0.526249e6 / 0.907200e6 -0.10819e5 / 0.777600e6 -0.50767e5 / 0.907200e6 -0.631e3 / 0.28800e5 0.91e2 / 0.7776e4; -0.526249e6 / 0.907200e6 0 0.1421209e7 / 0.2721600e7 0.16657e5 / 0.201600e6 -0.8467e4 / 0.453600e6 -0.33059e5 / 0.5443200e7; 0.10819e5 / 0.777600e6 -0.1421209e7 / 0.2721600e7 0 0.631187e6 / 0.1360800e7 0.400139e6 / 0.5443200e7 -0.8789e4 / 0.302400e6; 0.50767e5 / 0.907200e6 -0.16657e5 / 0.201600e6 -0.631187e6 / 0.1360800e7 0 0.496403e6 / 0.907200e6 -0.308533e6 / 0.5443200e7; 0.631e3 / 0.28800e5 0.8467e4 / 0.453600e6 -0.400139e6 / 0.5443200e7 -0.496403e6 / 0.907200e6 0 0.1805647e7 / 0.2721600e7; -0.91e2 / 0.7776e4 0.33059e5 / 0.5443200e7 0.8789e4 / 0.302400e6 0.308533e6 / 0.5443200e7 -0.1805647e7 / 0.2721600e7 0;];
    Q(1:6,1:6)=Q_U;
    Q(m-5:m,m-5:m)=rot90(  -Q_U ,2 );

    e_1=sparse(m,1);e_1(1)=1;
    e_m=sparse(m,1);e_m(m)=1;


    D1=H\(Q-1/2*(e_1*e_1')+1/2*(e_m*e_m')) ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    % Second derivative, 1st order accurate at first 6 boundary points
%     m2=1/12;m1=-16/12;m0=30/12;
%     M=m2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2))+m1*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1))+m0*diag(ones(m,1),0);
    %M=(1/12*diag(ones(m-2,1),2)-16/12*diag(ones(m-1,1),1)-16/12*diag(ones(m-1,1),-1)+1/12*diag(ones(m-2,1),-2)+30/12*diag(ones(m,1),0));
    M=-spdiags([-e 16*e -30*e 16*e -e], -2:2, m, m)/12;
    M_U=[0.2386127e7 / 0.2177280e7 -0.515449e6 / 0.453600e6 -0.10781e5 / 0.777600e6 0.61567e5 / 0.1360800e7 0.6817e4 / 0.403200e6 -0.1069e4 / 0.136080e6; -0.515449e6 / 0.453600e6 0.4756039e7 / 0.2177280e7 -0.1270009e7 / 0.1360800e7 -0.3751e4 / 0.28800e5 0.3067e4 / 0.680400e6 0.119459e6 / 0.10886400e8; -0.10781e5 / 0.777600e6 -0.1270009e7 / 0.1360800e7 0.111623e6 / 0.60480e5 -0.555587e6 / 0.680400e6 -0.551339e6 / 0.5443200e7 0.8789e4 / 0.453600e6; 0.61567e5 / 0.1360800e7 -0.3751e4 / 0.28800e5 -0.555587e6 / 0.680400e6 0.1025327e7 / 0.544320e6 -0.464003e6 / 0.453600e6 0.222133e6 / 0.5443200e7; 0.6817e4 / 0.403200e6 0.3067e4 / 0.680400e6 -0.551339e6 / 0.5443200e7 -0.464003e6 / 0.453600e6 0.5074159e7 / 0.2177280e7 -0.1784047e7 / 0.1360800e7; -0.1069e4 / 0.136080e6 0.119459e6 / 0.10886400e8 0.8789e4 / 0.453600e6 0.222133e6 / 0.5443200e7 -0.1784047e7 / 0.1360800e7 0.1812749e7 / 0.725760e6;];

    M(1:6,1:6)=M_U;

    M(m-5:m,m-5:m)=rot90(  M_U ,2 );
    M=M/h;

    S_U=[-0.11e2 / 0.6e1 3 -0.3e1 / 0.2e1 0.1e1 / 0.3e1;]/h;
    S_1=sparse(1,m);
    S_1(1:4)=S_U;
    S_m=sparse(1,m);

    S_m(m-3:m)=fliplr(-S_U);

    D2=H\(-M-e_1*S_1+e_m*S_m);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    % Third derivative, 1st order accurate at first 6 boundary points

    q3=-1/8;q2=1;q1=-13/8;
%     Q3=q3*(diag(ones(m-3,1),3)-diag(ones(m-3,1),-3))+q2*(diag(ones(m-2,1),2)-diag(ones(m-2,1),-2))+q1*(diag(ones(m-1,1),1)-diag(ones(m-1,1),-1));
    diags   = -3:3;
    stencil = [-q3,-q2,-q1,0,q1,q2,q3];
    Q3 = stripeMatrix(stencil, diags, m);

    %QQ3=(-1/8*diag(ones(m-3,1),3) + 1*diag(ones(m-2,1),2) - 13/8*diag(ones(m-1,1),1) +13/8*diag(ones(m-1,1),-1) -1*diag(ones(m-2,1),-2) + 1/8*diag(ones(m-3,1),-3));


    Q3_U = [0 -0.88471e5 / 0.67200e5 0.58139e5 / 0.33600e5 -0.1167e4 / 0.2800e4 -0.89e2 / 0.11200e5 0.7e1 / 0.640e3; 0.88471e5 / 0.67200e5 0 -0.43723e5 / 0.16800e5 0.46783e5 / 0.33600e5 -0.191e3 / 0.3200e4 -0.1567e4 / 0.33600e5; -0.58139e5 / 0.33600e5 0.43723e5 / 0.16800e5 0 -0.4049e4 / 0.2400e4 0.29083e5 / 0.33600e5 -0.71e2 / 0.1400e4; 0.1167e4 / 0.2800e4 -0.46783e5 / 0.33600e5 0.4049e4 / 0.2400e4 0 -0.8591e4 / 0.5600e4 0.10613e5 / 0.11200e5; 0.89e2 / 0.11200e5 0.191e3 / 0.3200e4 -0.29083e5 / 0.33600e5 0.8591e4 / 0.5600e4 0 -0.108271e6 / 0.67200e5; -0.7e1 / 0.640e3 0.1567e4 / 0.33600e5 0.71e2 / 0.1400e4 -0.10613e5 / 0.11200e5 0.108271e6 / 0.67200e5 0;];

    Q3(1:6,1:6)=Q3_U;
    Q3(m-5:m,m-5:m)=rot90(  -Q3_U ,2 );
    Q3=Q3/h^2;



    S2_U=[2 -5 4 -1;]/h^2;
    S2_1=sparse(1,m);
    S2_1(1:4)=S2_U;
    S2_m=sparse(1,m);
    S2_m(m-3:m)=fliplr(S2_U);



    D3=H\(Q3 - e_1*S2_1 + e_m*S2_m +1/2*(S_1'*S_1) -1/2*(S_m'*S_m) ) ;

    % Fourth derivative, 0th order accurate at first 6 boundary points (still
    % yield 4th order convergence if stable: for example u_tt=-u_xxxx

    m3=-1/6;m2=2;m1=-13/2;m0=28/3;
%     M4=m3*(diag(ones(m-3,1),3)+diag(ones(m-3,1),-3))+m2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2))+m1*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1))+m0*diag(ones(m,1),0);
    diags   = -3:3;
    left_stencil = [m3,m2,m1];
    stencil = [left_stencil,m0,fliplr(left_stencil)];
    M4 = stripeMatrix(stencil, diags, m);

    %M4=(-1/6*(diag(ones(m-3,1),3)+diag(ones(m-3,1),-3) ) + 2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2)) -13/2*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1)) + 28/3*diag(ones(m,1),0));

    M4_U=[0.4596181e7 / 0.1814400e7 -0.10307743e8 / 0.1814400e7 0.160961e6 / 0.43200e5 -0.535019e6 / 0.907200e6 0.109057e6 / 0.1814400e7 -0.29273e5 / 0.604800e6; -0.10307743e8 / 0.1814400e7 0.8368543e7 / 0.604800e6 -0.9558943e7 / 0.907200e6 0.2177057e7 / 0.907200e6 -0.11351e5 / 0.86400e5 0.204257e6 / 0.1814400e7; 0.160961e6 / 0.43200e5 -0.9558943e7 / 0.907200e6 0.4938581e7 / 0.453600e6 -0.786473e6 / 0.151200e6 0.1141057e7 / 0.907200e6 -0.120619e6 / 0.907200e6; -0.535019e6 / 0.907200e6 0.2177057e7 / 0.907200e6 -0.786473e6 / 0.151200e6 0.3146581e7 / 0.453600e6 -0.4614143e7 / 0.907200e6 0.24587e5 / 0.14400e5; 0.109057e6 / 0.1814400e7 -0.11351e5 / 0.86400e5 0.1141057e7 / 0.907200e6 -0.4614143e7 / 0.907200e6 0.185709e6 / 0.22400e5 -0.11293343e8 / 0.1814400e7; -0.29273e5 / 0.604800e6 0.204257e6 / 0.1814400e7 -0.120619e6 / 0.907200e6 0.24587e5 / 0.14400e5 -0.11293343e8 / 0.1814400e7 0.16787381e8 / 0.1814400e7;];

    M4(1:6,1:6)=M4_U;

    M4(m-5:m,m-5:m)=rot90(  M4_U ,2 );
    M4=M4/h^3;

    S3_U=[-1 3 -3 1;]/h^3;
    S3_1=sparse(1,m);
    S3_1(1:4)=S3_U;
    S3_m=sparse(1,m);
    S3_m(m-3:m)=fliplr(-S3_U);

    D4=H\(M4-e_1*S3_1+e_m*S3_m  + S_1'*S2_1-S_m'*S2_m);


    % L=h*(m-1);
    %
    % x1=linspace(0,L,m)';
    % x2=x1.^2/fac(2);
    % x3=x1.^3/fac(3);
    % x4=x1.^4/fac(4);
    % x5=x1.^5/fac(5);
    %
    % x0=x1.^0/fac(1);

    S_1  = S_1';
    S2_1 = S2_1';
    S3_1 = S3_1';
    S_m  = S_m';
    S2_m = S2_m';
    S3_m = S3_m';



end
