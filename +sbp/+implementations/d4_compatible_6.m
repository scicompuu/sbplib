function [H, HI, D1, D4, e_1, e_m, M4, Q, S2_1, S2_m,...
    S3_1, S3_m, S_1, S_m] = d4_compatible_6(m,h)
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
    
    BP = 8;
    if(m<2*BP)
        error(['Operator requires at least ' num2str(2*BP) ' grid points']);
    end

    H=speye(m,m);
    H_U=[0.7493827e7 / 0.25401600e8 0 0 0 0 0 0 0; 0 0.5534051e7 / 0.3628800e7 0 0 0 0 0 0; 0 0 0.104561e6 / 0.403200e6 0 0 0 0 0; 0 0 0 0.260503e6 / 0.145152e6 0 0 0 0; 0 0 0 0 0.43237e5 / 0.103680e6 0 0 0; 0 0 0 0 0 0.514081e6 / 0.403200e6 0 0; 0 0 0 0 0 0 0.3356179e7 / 0.3628800e7 0; 0 0 0 0 0 0 0 0.25631027e8 / 0.25401600e8;];

    H(1:8,1:8)=H_U;
    H(m-7:m,m-7:m)=rot90(H_U,2);
    H=H*h;
    HI=inv(H);


    % First derivative SBP operator, 3rd order accurate at first 8 boundary points

    q3=1/60;q2=-3/20;q1=3/4;
%     Q=q3*(diag(ones(m-3,1),3) - diag(ones(m-3,1),-3))+q2*(diag(ones(m-2,1),2) - diag(ones(m-2,1),-2))+q1*(diag(ones(m-1,1),1)-diag(ones(m-1,1),-1));
    stencil = [-q3,-q2,-q1,0,q1,q2,q3];
    d = (length(stencil)-1)/2;
    diags = -d:d;
    Q = stripeMatrix(stencil, diags, m);

    Q_U = [0 0.26431903e8 / 0.43545600e8 0.39791489e8 / 0.152409600e9 -0.12747751e8 / 0.16934400e8 0.76099447e8 / 0.152409600e9 -0.1397443e7 / 0.12192768e8 0 0; -0.26431903e8 / 0.43545600e8 0 -0.13847476213e11 / 0.19559232000e11 0.35844843977e11 / 0.11735539200e11 -0.63413503537e11 / 0.23471078400e11 0.4764412871e10 / 0.3911846400e10 -0.1668252557e10 / 0.5867769600e10 0.842644697e9 / 0.29338848000e11; -0.39791489e8 / 0.152409600e9 0.13847476213e11 / 0.19559232000e11 0 -0.73834802771e11 / 0.23471078400e11 0.1802732209e10 / 0.325987200e9 -0.65514173e8 / 0.16299360e8 0.79341409141e11 / 0.58677696000e11 -0.1282384321e10 / 0.7823692800e10; 0.12747751e8 / 0.16934400e8 -0.35844843977e11 / 0.11735539200e11 0.73834802771e11 / 0.23471078400e11 0 -0.5274106087e10 / 0.1173553920e10 0.33743985841e11 / 0.5867769600e10 -0.6482602549e10 / 0.2607897600e10 0.1506017269e10 / 0.3911846400e10; -0.76099447e8 / 0.152409600e9 0.63413503537e11 / 0.23471078400e11 -0.1802732209e10 / 0.325987200e9 0.5274106087e10 / 0.1173553920e10 0 -0.7165829063e10 / 0.2607897600e10 0.23903110999e11 / 0.11735539200e11 -0.5346675911e10 / 0.11735539200e11; 0.1397443e7 / 0.12192768e8 -0.4764412871e10 / 0.3911846400e10 0.65514173e8 / 0.16299360e8 -0.33743985841e11 / 0.5867769600e10 0.7165829063e10 / 0.2607897600e10 0 -0.1060918223e10 / 0.11735539200e11 0.628353989e9 / 0.3911846400e10; 0 0.1668252557e10 / 0.5867769600e10 -0.79341409141e11 / 0.58677696000e11 0.6482602549e10 / 0.2607897600e10 -0.23903110999e11 / 0.11735539200e11 0.1060918223e10 / 0.11735539200e11 0 0.25889988599e11 / 0.39118464000e11; 0 -0.842644697e9 / 0.29338848000e11 0.1282384321e10 / 0.7823692800e10 -0.1506017269e10 / 0.3911846400e10 0.5346675911e10 / 0.11735539200e11 -0.628353989e9 / 0.3911846400e10 -0.25889988599e11 / 0.39118464000e11 0;];

    Q(1:8,1:8)=Q_U;
    Q(m-7:m,m-7:m)=rot90(  -Q_U ,2 );

    e_1=sparse(m,1);e_1(1)=1;
    e_m=sparse(m,1);e_m(m)=1;


    D1=H\(Q-1/2*(e_1*e_1')+1/2*(e_m*e_m')) ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    % Second derivative, 1st order accurate at first 6 boundary points
    % m3=-1/90;m2=3/20;m1=-3/2;m0=49/18;
    %
    % M=m3*(diag(ones(m-3,1),3)+diag(ones(m-3,1),-3))+m2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2))+m1*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1))+m0*diag(ones(m,1),0);
    % M_U=[0.4347276223e10 / 0.3736212480e10 -0.1534657609e10 / 0.1210809600e10 0.68879e5 / 0.3057600e7 0.1092927401e10 / 0.13076743680e11 0.18145423e8 / 0.968647680e9 -0.1143817e7 / 0.60540480e8 -0.355447739e9 / 0.65383718400e11 0.56081e5 / 0.16473600e8; -0.1534657609e10 / 0.1210809600e10 0.42416226217e11 / 0.18681062400e11 -0.228654119e9 / 0.345945600e9 -0.12245627e8 / 0.34594560e8 -0.2995295e7 / 0.46702656e8 0.52836503e8 / 0.691891200e9 0.119351e6 / 0.12812800e8 -0.634102039e9 / 0.65383718400e11; 0.68879e5 / 0.3057600e7 -0.228654119e9 / 0.345945600e9 0.5399287e7 / 0.4193280e7 -0.24739409e8 / 0.34594560e8 0.7878737e7 / 0.69189120e8 -0.1917829e7 / 0.31449600e8 0.39727e5 / 0.3660800e7 0.10259e5 / 0.4656960e7; 0.1092927401e10 / 0.13076743680e11 -0.12245627e8 / 0.34594560e8 -0.24739409e8 / 0.34594560e8 0.7780367599e10 / 0.3736212480e10 -0.70085363e8 / 0.69189120e8 -0.500209e6 / 0.6289920e7 -0.311543e6 / 0.17962560e8 0.278191e6 / 0.21525504e8; 0.18145423e8 / 0.968647680e9 -0.2995295e7 / 0.46702656e8 0.7878737e7 / 0.69189120e8 -0.70085363e8 / 0.69189120e8 0.7116321131e10 / 0.3736212480e10 -0.545081e6 / 0.532224e6 0.811631e6 / 0.11531520e8 -0.84101639e8 / 0.13076743680e11; -0.1143817e7 / 0.60540480e8 0.52836503e8 / 0.691891200e9 -0.1917829e7 / 0.31449600e8 -0.500209e6 / 0.6289920e7 -0.545081e6 / 0.532224e6 0.324760747e9 / 0.138378240e9 -0.65995697e8 / 0.49420800e8 0.1469203e7 / 0.13759200e8; -0.355447739e9 / 0.65383718400e11 0.119351e6 / 0.12812800e8 0.39727e5 / 0.3660800e7 -0.311543e6 / 0.17962560e8 0.811631e6 / 0.11531520e8 -0.65995697e8 / 0.49420800e8 0.48284442317e11 / 0.18681062400e11 -0.1762877569e10 / 0.1210809600e10; 0.56081e5 / 0.16473600e8 -0.634102039e9 / 0.65383718400e11 0.10259e5 / 0.4656960e7 0.278191e6 / 0.21525504e8 -0.84101639e8 / 0.13076743680e11 0.1469203e7 / 0.13759200e8 -0.1762877569e10 / 0.1210809600e10 0.10117212851e11 / 0.3736212480e10;];
    %
    % M(1:8,1:8)=M_U;
    %
    % M(m-7:m,m-7:m)=flipud( fliplr( M_U ) );
    % M=M/h;

    S_U=[-0.12700800e8 / 0.7493827e7 0.185023321e9 / 0.89925924e8 0.39791489e8 / 0.44962962e8 -0.38243253e8 / 0.14987654e8 0.76099447e8 / 0.44962962e8 -0.34936075e8 / 0.89925924e8;]/h;
    S_1=sparse(1,m);
    S_1(1:6)=S_U;
    S_m=sparse(1,m);

    S_m(m-5:m)=fliplr(-S_U);

    %D2=H\(-M-e_1*S_1+e_m*S_m);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    % Third derivative, 1st order accurate at first 6 boundary points

    % q4=7/240;q3=-3/10;q2=169/120;q1=-61/30;
    % Q3=q4*(diag(ones(m-4,1),4)-diag(ones(m-4,1),-4))+q3*(diag(ones(m-3,1),3)-diag(ones(m-3,1),-3))+q2*(diag(ones(m-2,1),2)-diag(ones(m-2,1),-2))+q1*(diag(ones(m-1,1),1)-diag(ones(m-1,1),-1));
    %
    % %QQ3=(-1/8*diag(ones(m-3,1),3) + 1*diag(ones(m-2,1),2) - 13/8*diag(ones(m-1,1),1) +13/8*diag(ones(m-1,1),-1) -1*diag(ones(m-2,1),-2) + 1/8*diag(ones(m-3,1),-3));
    %
    %
    % Q3_U = [0 -0.10882810591e11 / 0.5811886080e10 0.398713069e9 / 0.132088320e9 -0.1746657571e10 / 0.1162377216e10 0.56050639e8 / 0.145297152e9 -0.11473393e8 / 0.1162377216e10 -0.38062741e8 / 0.1452971520e10 0.30473e5 / 0.4392960e7; 0.10882810591e11 / 0.5811886080e10 0 -0.3720544343e10 / 0.830269440e9 0.767707019e9 / 0.207567360e9 -0.1047978301e10 / 0.830269440e9 0.1240729e7 / 0.14826240e8 0.6807397e7 / 0.55351296e8 -0.50022767e8 / 0.1452971520e10; -0.398713069e9 / 0.132088320e9 0.3720544343e10 / 0.830269440e9 0 -0.2870078009e10 / 0.830269440e9 0.74962049e8 / 0.29652480e8 -0.12944857e8 / 0.30750720e8 -0.17846623e8 / 0.103783680e9 0.68707591e8 / 0.1162377216e10; 0.1746657571e10 / 0.1162377216e10 -0.767707019e9 / 0.207567360e9 0.2870078009e10 / 0.830269440e9 0 -0.727867087e9 / 0.276756480e9 0.327603877e9 / 0.207567360e9 -0.175223717e9 / 0.830269440e9 0.1353613e7 / 0.726485760e9; -0.56050639e8 / 0.145297152e9 0.1047978301e10 / 0.830269440e9 -0.74962049e8 / 0.29652480e8 0.727867087e9 / 0.276756480e9 0 -0.1804641793e10 / 0.830269440e9 0.311038417e9 / 0.207567360e9 -0.1932566239e10 / 0.5811886080e10; 0.11473393e8 / 0.1162377216e10 -0.1240729e7 / 0.14826240e8 0.12944857e8 / 0.30750720e8 -0.327603877e9 / 0.207567360e9 0.1804641793e10 / 0.830269440e9 0 -0.1760949511e10 / 0.830269440e9 0.2105883973e10 / 0.1452971520e10; 0.38062741e8 / 0.1452971520e10 -0.6807397e7 / 0.55351296e8 0.17846623e8 / 0.103783680e9 0.175223717e9 / 0.830269440e9 -0.311038417e9 / 0.207567360e9 0.1760949511e10 / 0.830269440e9 0 -0.1081094773e10 / 0.528353280e9; -0.30473e5 / 0.4392960e7 0.50022767e8 / 0.1452971520e10 -0.68707591e8 / 0.1162377216e10 -0.1353613e7 / 0.726485760e9 0.1932566239e10 / 0.5811886080e10 -0.2105883973e10 / 0.1452971520e10 0.1081094773e10 / 0.528353280e9 0;];
    %
    % Q3(1:8,1:8)=Q3_U;
    % Q3(m-7:m,m-7:m)=flipud( fliplr( -Q3_U ) );
    % Q3=Q3/h^2;
    %
    %
    %
     S2_U=[0.35e2 / 0.12e2 -0.26e2 / 0.3e1 0.19e2 / 0.2e1 -0.14e2 / 0.3e1 0.11e2 / 0.12e2;]/h^2;
     S2_1=sparse(1,m);
     S2_1(1:5)=S2_U;
     S2_m=sparse(1,m);
     S2_m(m-4:m)=fliplr(S2_U);
    %
    %
    %
    % D3=H\(Q3 - e_1*S2_1 + e_m*S2_m +1/2*(S_1'*S_1) -1/2*(S_m'*S_m) ) ;

    % Fourth derivative, 0th order accurate at first 6 boundary points (still
    % yield 4th order convergence if stable: for example u_tt=-u_xxxx

    m4=7/240;m3=-2/5;m2=169/60;m1=-122/15;m0=91/8;
%     M4=m4*(diag(ones(m-4,1),4)+diag(ones(m-4,1),-4))+m3*(diag(ones(m-3,1),3)+diag(ones(m-3,1),-3))+m2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2))+m1*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1))+m0*diag(ones(m,1),0);
    stencil = [m4,m3,m2,m1,m0,m1,m2,m3,m4];
    d = (length(stencil)-1)/2;
    diags = -d:d;
    M4 = stripeMatrix(stencil, diags, m);


    M4_U=[0.600485868980522851e18 / 0.274825314114120000e18 -0.1421010223681841e16 / 0.348984525859200e15 0.38908412970187e14 / 0.1586293299360000e16 0.10224077451922837e17 / 0.2243471951952000e16 -0.7577302712815639e16 / 0.1744922629296000e16 0.138091642084013e15 / 0.59351109840000e14 -0.3775041725375197e16 / 0.4486943903904000e16 0.9907210230881393e16 / 0.61072292025360000e17; -0.1421010223681841e16 / 0.348984525859200e15 0.3985852497808703e16 / 0.407903991264000e15 -0.90048788923861e14 / 0.15579666333000e14 -0.4312795866499e13 / 0.997098645312e12 0.4414634708891947e16 / 0.448694390390400e15 -0.886174803100459e15 / 0.99709864531200e14 0.4333e4 / 0.1000e4 -0.13800578064893047e17 / 0.15704303663664000e17; 0.38908412970187e14 / 0.1586293299360000e16 -0.90048788923861e14 / 0.15579666333000e14 0.2071682582321887e16 / 0.113306664240000e15 -0.769471337294003e15 / 0.41545776888000e14 0.112191585452033e15 / 0.166183107552000e15 0.7204491902193671e16 / 0.623186653320000e15 -0.24847093554379e14 / 0.3115933266600e13 0.943854037768721e15 / 0.545288321655000e15; 0.10224077451922837e17 / 0.2243471951952000e16 -0.4312795866499e13 / 0.997098645312e12 -0.769471337294003e15 / 0.41545776888000e14 0.3086874339649421e16 / 0.81580798252800e14 -0.396009005312111e15 / 0.16618310755200e14 0.348854811893087e15 / 0.249274661328000e15 0.895954627955053e15 / 0.224347195195200e15 -0.184881685054543e15 / 0.166183107552000e15; -0.7577302712815639e16 / 0.1744922629296000e16 0.4414634708891947e16 / 0.448694390390400e15 0.112191585452033e15 / 0.166183107552000e15 -0.396009005312111e15 / 0.16618310755200e14 0.3774861828677557e16 / 0.112173597597600e15 -0.5693689108983593e16 / 0.249274661328000e15 0.803944126167107e15 / 0.99709864531200e14 -0.19547569411550791e17 / 0.15704303663664000e17; 0.138091642084013e15 / 0.59351109840000e14 -0.886174803100459e15 / 0.99709864531200e14 0.7204491902193671e16 / 0.623186653320000e15 0.348854811893087e15 / 0.249274661328000e15 -0.5693689108983593e16 / 0.249274661328000e15 0.73965842628398389e17 / 0.2492746613280000e16 -0.2184472662036043e16 / 0.124637330664000e15 0.46667e5 / 0.10000e5; -0.3775041725375197e16 / 0.4486943903904000e16 0.4333e4 / 0.1000e4 -0.24847093554379e14 / 0.3115933266600e13 0.895954627955053e15 / 0.224347195195200e15 0.803944126167107e15 / 0.99709864531200e14 -0.2184472662036043e16 / 0.124637330664000e15 0.37593640125444199e17 / 0.2243471951952000e16 -0.37e2 / 0.4e1; 0.9907210230881393e16 / 0.61072292025360000e17 -0.13800578064893047e17 / 0.15704303663664000e17 0.943854037768721e15 / 0.545288321655000e15 -0.184881685054543e15 / 0.166183107552000e15 -0.19547569411550791e17 / 0.15704303663664000e17 0.46667e5 / 0.10000e5 -0.37e2 / 0.4e1 0.12766926490502478779e20 / 0.1099301256456480000e19;];

    M4(1:8,1:8)=M4_U;

    M4(m-7:m,m-7:m)=rot90(  M4_U ,2 );
    M4=M4/h^3;

    S3_U=[-0.5e1 / 0.2e1 9 -12 7 -0.3e1 / 0.2e1;]/h^3;
    S3_1=sparse(1,m);
    S3_1(1:5)=S3_U;
    S3_m=sparse(1,m);
    S3_m(m-4:m)=fliplr(-S3_U);

    D4=H\(M4-e_1*S3_1+e_m*S3_m  + S_1'*S2_1-S_m'*S2_m);

    S_1 = S_1';
    S_m = S_m';
    S2_1 = S2_1';
    S2_m = S2_m';
    S3_1 = S3_1';
    S3_m = S3_m';

    % L=h*(m-1);

    % x1=linspace(0,L,m)';
    % x2=x1.^2/fac(2);
    % x3=x1.^3/fac(3);
    % x4=x1.^4/fac(4);
    % x5=x1.^5/fac(5);

    % x0=x1.^0/fac(1);

end
