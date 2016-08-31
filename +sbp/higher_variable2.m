% Returns D2 as a function handle
function [H, HI, D1, D2, D3, D4, e_1, e_m, M4, Q, S2_1, S2_m, S3_1, S3_m, S_1, S_m] = higher_variable2(m,h)
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

    H=diag(ones(m,1),0);H(1,1)=1/2;H(m,m)=1/2;


    H=H*h;
    HI=inv(H);


    % First derivative SBP operator, 1st order accurate at first 6 boundary points

    q1=1/2;
    Q=q1*(diag(ones(m-1,1),1)-diag(ones(m-1,1),-1));

    %Q=(-1/12*diag(ones(m-2,1),2)+8/12*diag(ones(m-1,1),1)-8/12*diag(ones(m-1,1),-1)+1/12*diag(ones(m-2,1),-2));


    e_1=zeros(m,1);e_1(1)=1;
    e_m=zeros(m,1);e_m(m)=1;


    D1=HI*(Q-1/2*e_1*e_1'+1/2*e_m*e_m') ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    % Second derivative, 1st order accurate at first boundary points

    %% below for constant coefficients
    % m1=-1;m0=2;
    % M=m1*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1))+m0*diag(ones(m,1),0);M(1,1)=1;M(m,m)=1;
    % M=M/h;
    %D2=HI*(-M-e_1*S_1+e_m*S_m);

    %% Below for variable coefficients
    %% Require a vector c with the koeffients

    S_U=[-3/2 2 -1/2]/h;
    S_1=zeros(1,m);
    S_1(1:3)=S_U;
    S_m=zeros(1,m);
    S_m(m-2:m)=fliplr(-S_U);

    S_1 = S_1';
    S_m = S_m';

    M=sparse(m,m);
    e_1 = sparse(e_1);
    e_m = sparse(e_m);
    S_1 = sparse(S_1);
    S_m = sparse(S_m);

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

        D2=HI*(-M-c(1)*e_1*S_1'+c(m)*e_m*S_m');
    end
    D2 = @D2_fun;





    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    % Third derivative, 1st order accurate at first 6 boundary points

    q2=1/2;q1=-1;
    Q3=q2*(diag(ones(m-2,1),2)-diag(ones(m-2,1),-2))+q1*(diag(ones(m-1,1),1)-diag(ones(m-1,1),-1));

    %QQ3=(-1/8*diag(ones(m-3,1),3) + 1*diag(ones(m-2,1),2) - 13/8*diag(ones(m-1,1),1) +13/8*diag(ones(m-1,1),-1) -1*diag(ones(m-2,1),-2) + 1/8*diag(ones(m-3,1),-3));


    Q3_U = [0 -0.13e2 / 0.16e2 0.7e1 / 0.8e1 -0.1e1 / 0.16e2; 0.13e2 / 0.16e2 0 -0.23e2 / 0.16e2 0.5e1 / 0.8e1; -0.7e1 / 0.8e1 0.23e2 / 0.16e2 0 -0.17e2 / 0.16e2; 0.1e1 / 0.16e2 -0.5e1 / 0.8e1 0.17e2 / 0.16e2 0;];
    Q3(1:4,1:4)=Q3_U;
    Q3(m-3:m,m-3:m)=flipud( fliplr( -Q3_U ) );
    Q3=Q3/h^2;



    S2_U=[1 -2 1;]/h^2;
    S2_1=zeros(1,m);
    S2_1(1:3)=S2_U;
    S2_m=zeros(1,m);
    S2_m(m-2:m)=fliplr(S2_U);
    S2_1 = S2_1';
    S2_m = S2_m';



    D3=HI*(Q3 - e_1*S2_1' + e_m*S2_m' +1/2*S_1*S_1' -1/2*S_m*S_m' ) ;

    % Fourth derivative, 0th order accurate at first 6 boundary points (still
    % yield 4th order convergence if stable: for example u_tt=-u_xxxx

    m2=1;m1=-4;m0=6;
    M4=m2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2))+m1*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1))+m0*diag(ones(m,1),0);

    %M4=(-1/6*(diag(ones(m-3,1),3)+diag(ones(m-3,1),-3) ) + 2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2)) -13/2*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1)) + 28/3*diag(ones(m,1),0));

    M4_U=[0.13e2 / 0.10e2 -0.12e2 / 0.5e1 0.9e1 / 0.10e2 0.1e1 / 0.5e1; -0.12e2 / 0.5e1 0.26e2 / 0.5e1 -0.16e2 / 0.5e1 0.2e1 / 0.5e1; 0.9e1 / 0.10e2 -0.16e2 / 0.5e1 0.47e2 / 0.10e2 -0.17e2 / 0.5e1; 0.1e1 / 0.5e1 0.2e1 / 0.5e1 -0.17e2 / 0.5e1 0.29e2 / 0.5e1;];


    M4(1:4,1:4)=M4_U;

    M4(m-3:m,m-3:m)=flipud( fliplr( M4_U ) );
    M4=M4/h^3;

    S3_U=[-1 3 -3 1;]/h^3;
    S3_1=zeros(1,m);
    S3_1(1:4)=S3_U;
    S3_m=zeros(1,m);
    S3_m(m-3:m)=fliplr(-S3_U);
    S3_1 = S3_1';
    S3_m = S3_m';

    D4=HI*(M4-e_1*S3_1'+e_m*S3_m'  + S_1*S2_1'-S_m*S2_m');
end
