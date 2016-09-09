function [H, HI, D1, D4, e_1, e_m, M4, Q, S2_1, S2_m,...
    S3_1, S3_m, S_1, S_m] = d4_compatible_4(m,h)
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

    H=diag(ones(m,1),0);
    H_U=[0.3e1 / 0.11e2 0 0 0 0 0; 0 0.2125516311e10 / 0.1311004640e10 0 0 0 0; 0 0 0.278735189e9 / 0.1966506960e10 0 0 0; 0 0 0 0.285925927e9 / 0.163875580e9 0 0; 0 0 0 0 0.1284335339e10 / 0.1966506960e10 0; 0 0 0 0 0 0.4194024163e10 / 0.3933013920e10;];
    H(1:6,1:6)=H_U;
    H(m-5:m,m-5:m)=fliplr(flipud(H_U));
    H=H*h;
    HI=inv(H);


    % First derivative SBP operator, 1st order accurate at first 6 boundary points

    q2=-1/12;q1=8/12;
    Q=q2*(diag(ones(m-2,1),2) - diag(ones(m-2,1),-2))+q1*(diag(ones(m-1,1),1)-diag(ones(m-1,1),-1));

    %Q=(-1/12*diag(ones(m-2,1),2)+8/12*diag(ones(m-1,1),1)-8/12*diag(ones(m-1,1),-1)+1/12*diag(ones(m-2,1),-2));

    Q_U = [0 0.9e1 / 0.11e2 -0.9e1 / 0.22e2 0.1e1 / 0.11e2 0 0; -0.9e1 / 0.11e2 0 0.2595224893e10 / 0.2622009280e10 -0.151435707e9 / 0.327751160e9 0.1112665611e10 / 0.2622009280e10 -0.1290899e7 / 0.9639740e7; 0.9e1 / 0.22e2 -0.2595224893e10 / 0.2622009280e10 0 0.1468436423e10 / 0.983253480e9 -0.1194603401e10 / 0.983253480e9 0.72033031e8 / 0.238364480e9; -0.1e1 / 0.11e2 0.151435707e9 / 0.327751160e9 -0.1468436423e10 / 0.983253480e9 0 0.439819541e9 / 0.327751160e9 -0.215942641e9 / 0.983253480e9; 0 -0.1112665611e10 / 0.2622009280e10 0.1194603401e10 / 0.983253480e9 -0.439819541e9 / 0.327751160e9 0 0.1664113643e10 / 0.2622009280e10; 0 0.1290899e7 / 0.9639740e7 -0.72033031e8 / 0.238364480e9 0.215942641e9 / 0.983253480e9 -0.1664113643e10 / 0.2622009280e10 0;];
    Q(1:6,1:6)=Q_U;
    Q(m-5:m,m-5:m)=flipud( fliplr( -Q_U ) );

    e_1=zeros(m,1);e_1(1)=1;
    e_m=zeros(m,1);e_m(m)=1;


    D1=HI*(Q-1/2*e_1*e_1'+1/2*e_m*e_m') ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    % % Second derivative, 1st order accurate at first 6 boundary points
    % m2=1/12;m1=-16/12;m0=30/12;
    % M=m2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2))+m1*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1))+m0*diag(ones(m,1),0);
    % %M=(1/12*diag(ones(m-2,1),2)-16/12*diag(ones(m-1,1),1)-16/12*diag(ones(m-1,1),-1)+1/12*diag(ones(m-2,1),-2)+30/12*diag(ones(m,1),0));
    % M_U=[0.2386127e7 / 0.2177280e7 -0.515449e6 / 0.453600e6 -0.10781e5 / 0.777600e6 0.61567e5 / 0.1360800e7 0.6817e4 / 0.403200e6 -0.1069e4 / 0.136080e6; -0.515449e6 / 0.453600e6 0.4756039e7 / 0.2177280e7 -0.1270009e7 / 0.1360800e7 -0.3751e4 / 0.28800e5 0.3067e4 / 0.680400e6 0.119459e6 / 0.10886400e8; -0.10781e5 / 0.777600e6 -0.1270009e7 / 0.1360800e7 0.111623e6 / 0.60480e5 -0.555587e6 / 0.680400e6 -0.551339e6 / 0.5443200e7 0.8789e4 / 0.453600e6; 0.61567e5 / 0.1360800e7 -0.3751e4 / 0.28800e5 -0.555587e6 / 0.680400e6 0.1025327e7 / 0.544320e6 -0.464003e6 / 0.453600e6 0.222133e6 / 0.5443200e7; 0.6817e4 / 0.403200e6 0.3067e4 / 0.680400e6 -0.551339e6 / 0.5443200e7 -0.464003e6 / 0.453600e6 0.5074159e7 / 0.2177280e7 -0.1784047e7 / 0.1360800e7; -0.1069e4 / 0.136080e6 0.119459e6 / 0.10886400e8 0.8789e4 / 0.453600e6 0.222133e6 / 0.5443200e7 -0.1784047e7 / 0.1360800e7 0.1812749e7 / 0.725760e6;];
    %
    % M(1:6,1:6)=M_U;
    %
    % M(m-5:m,m-5:m)=flipud( fliplr( M_U ) );
    % M=M/h;
    %
     S_U=[-0.11e2 / 0.6e1 3 -0.3e1 / 0.2e1 0.1e1 / 0.3e1;]/h;
     S_1=zeros(1,m);
     S_1(1:4)=S_U;
     S_m=zeros(1,m);

     S_m(m-3:m)=fliplr(-S_U);

    % D2=HI*(-M-e_1*S_1+e_m*S_m);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    % Third derivative, 1st order accurate at first 6 boundary points

    % q3=-1/8;q2=1;q1=-13/8;
    % Q3=q3*(diag(ones(m-3,1),3)-diag(ones(m-3,1),-3))+q2*(diag(ones(m-2,1),2)-diag(ones(m-2,1),-2))+q1*(diag(ones(m-1,1),1)-diag(ones(m-1,1),-1));
    %
    % %QQ3=(-1/8*diag(ones(m-3,1),3) + 1*diag(ones(m-2,1),2) - 13/8*diag(ones(m-1,1),1) +13/8*diag(ones(m-1,1),-1) -1*diag(ones(m-2,1),-2) + 1/8*diag(ones(m-3,1),-3));
    %
    %
    % Q3_U = [0 -0.88471e5 / 0.67200e5 0.58139e5 / 0.33600e5 -0.1167e4 / 0.2800e4 -0.89e2 / 0.11200e5 0.7e1 / 0.640e3; 0.88471e5 / 0.67200e5 0 -0.43723e5 / 0.16800e5 0.46783e5 / 0.33600e5 -0.191e3 / 0.3200e4 -0.1567e4 / 0.33600e5; -0.58139e5 / 0.33600e5 0.43723e5 / 0.16800e5 0 -0.4049e4 / 0.2400e4 0.29083e5 / 0.33600e5 -0.71e2 / 0.1400e4; 0.1167e4 / 0.2800e4 -0.46783e5 / 0.33600e5 0.4049e4 / 0.2400e4 0 -0.8591e4 / 0.5600e4 0.10613e5 / 0.11200e5; 0.89e2 / 0.11200e5 0.191e3 / 0.3200e4 -0.29083e5 / 0.33600e5 0.8591e4 / 0.5600e4 0 -0.108271e6 / 0.67200e5; -0.7e1 / 0.640e3 0.1567e4 / 0.33600e5 0.71e2 / 0.1400e4 -0.10613e5 / 0.11200e5 0.108271e6 / 0.67200e5 0;];
    %
    % Q3(1:6,1:6)=Q3_U;
    % Q3(m-5:m,m-5:m)=flipud( fliplr( -Q3_U ) );
    % Q3=Q3/h^2;



    S2_U=[2 -5 4 -1;]/h^2;
    S2_1=zeros(1,m);
    S2_1(1:4)=S2_U;
    S2_m=zeros(1,m);
    S2_m(m-3:m)=fliplr(S2_U);



    %D3=HI*(Q3 - e_1*S2_1 + e_m*S2_m +1/2*S_1'*S_1 -1/2*S_m'*S_m ) ;

    % Fourth derivative, 0th order accurate at first 6 boundary points (still
    % yield 4th order convergence if stable: for example u_tt=-u_xxxx

    m3=-1/6;m2=2;m1=-13/2;m0=28/3;
    M4=m3*(diag(ones(m-3,1),3)+diag(ones(m-3,1),-3))+m2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2))+m1*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1))+m0*diag(ones(m,1),0);

    %M4=(-1/6*(diag(ones(m-3,1),3)+diag(ones(m-3,1),-3) ) + 2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2)) -13/2*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1)) + 28/3*diag(ones(m,1),0));

    M4_U=[0.227176919517319e15 / 0.94899692875680e14 -0.15262605263734e14 / 0.2965615402365e13 0.20205404771243e14 / 0.6778549491120e13 -0.3998303664097e13 / 0.23724923218920e14 0.1088305091927e13 / 0.94899692875680e14 -0.1686077077693e13 / 0.23724923218920e14; -0.15262605263734e14 / 0.2965615402365e13 0.280494781164181e15 / 0.23724923218920e14 -0.46417445546261e14 / 0.5931230804730e13 0.1705307929429e13 / 0.1694637372780e13 -0.553547394061e12 / 0.5931230804730e13 0.5615721694973e13 / 0.23724923218920e14; 0.20205404771243e14 / 0.6778549491120e13 -0.46417445546261e14 / 0.5931230804730e13 0.4135802350237e13 / 0.551742400440e12 -0.4140981465247e13 / 0.1078405600860e13 0.75538453067437e14 / 0.47449846437840e14 -0.4778134936391e13 / 0.11862461609460e14; -0.3998303664097e13 / 0.23724923218920e14 0.1705307929429e13 / 0.1694637372780e13 -0.4140981465247e13 / 0.1078405600860e13 0.20760974175677e14 / 0.2965615402365e13 -0.138330689701889e15 / 0.23724923218920e14 0.23711317526909e14 / 0.11862461609460e14; 0.1088305091927e13 / 0.94899692875680e14 -0.553547394061e12 / 0.5931230804730e13 0.75538453067437e14 / 0.47449846437840e14 -0.138330689701889e15 / 0.23724923218920e14 0.120223780251937e15 / 0.13557098982240e14 -0.151383731537477e15 / 0.23724923218920e14; -0.1686077077693e13 / 0.23724923218920e14 0.5615721694973e13 / 0.23724923218920e14 -0.4778134936391e13 / 0.11862461609460e14 0.23711317526909e14 / 0.11862461609460e14 -0.151383731537477e15 / 0.23724923218920e14 0.220304030094121e15 / 0.23724923218920e14;];

    M4(1:6,1:6)=M4_U;

    M4(m-5:m,m-5:m)=flipud( fliplr( M4_U ) );
    M4=M4/h^3;

    S3_U=[-1 3 -3 1;]/h^3;
    S3_1=zeros(1,m);
    S3_1(1:4)=S3_U;
    S3_m=zeros(1,m);
    S3_m(m-3:m)=fliplr(-S3_U);

    D4=HI*(M4-e_1*S3_1+e_m*S3_m  + S_1'*S2_1-S_m'*S2_m);

    S_1 = S_1';
    S_m = S_m';
    S2_1 = S2_1';
    S2_m = S2_m';
    S3_1 = S3_1';
    S3_m = S3_m';

    % L=h*(m-1);
    %
    % x1=linspace(0,L,m)';
    % x2=x1.^2/fac(2);
    % x3=x1.^3/fac(3);
    % x4=x1.^4/fac(4);
    % x5=x1.^5/fac(5);
    %
    % x0=x1.^0/fac(1);

end
