function [H, HI, D1, D2, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = d4_variable_6_2(m,h)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% 6:te ordn. SBP Finita differens         %%%
    %%% operatorer med diagonal norm            %%%
    %%% Extension to variable koeff             %%%
    %%%                                         %%%
    %%% H           (Normen)                    %%%
    %%% D1=H^(-1)Q  (approx f?rsta derivatan)   %%%
    %%% D2          (approx andra derivatan)    %%%
    %%% D2=HI*(R+C*D*S                          %%%
    %%%                                         %%%
    %%% R=-D1'*H*C*D1-RR                        %%%
    %%%                                         %%%
    %%% RR ?r dissipation)                      %%%
    %%% Dissipationen uppbyggd av D4:           %%%
    %%% DI=D4*B*H*D4                            %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % H?r med 6 RP ist?llet f?r 8 f?r D4 operatorn, dock samma randderivator
    % Denna ?r noggrannare, och har 2a ordningens randdslutning och b?r ge 6te
    % ordningens konvergens. Hade dock ingen fri parameter att optimera


    H = diag(ones(m,1),0);
    H(1:6,1:6) = [
        0.181e3/0.576e3 0 0 0 0 0;
        0 0.1343e4/0.960e3 0 0 0 0;
        0 0 0.293e3/0.480e3 0 0 0;
        0 0 0 0.1811e4/0.1440e4 0 0;
        0 0 0 0 0.289e3/0.320e3 0;
        0 0 0 0 0 0.65e2/0.64e2;
    ];
    H(m-5:m,m-5:m) = fliplr(flipud(H(1:6,1:6)));


    e_1 = zeros(m,1);e_1(1) = 1;
    e_m = zeros(m,1);e_m(m) = 1;

    S_U = [-0.137e3/0.60e2 5 -5 0.10e2/0.3e1 -0.5e1/0.4e1 0.1e1/0.5e1;]/h;
    S_1 = zeros(1,m);
    S_1(1:6) = S_U;
    S_m = zeros(1,m);
    S_m(m-5:m) = fliplr(-S_U);

    S2_U = [0.15e2/0.4e1 -0.77e2/0.6e1 0.107e3/0.6e1 -13 0.61e2/0.12e2 -0.5e1/0.6e1;]/h^2;
    S2_1 = zeros(1,m);
    S2_1(1:6) = S2_U;
    S2_m = zeros(1,m);
    S2_m(m-5:m) = fliplr(S2_U);

    S3_U = [-0.17e2/0.4e1 0.71e2/0.4e1 -0.59e2/0.2e1 0.49e2/0.2e1 -0.41e2/0.4e1 0.7e1/0.4e1;]/h^3;
    S3_1 = zeros(1,m);
    S3_1(1:6) = S3_U;
    S3_m = zeros(1,m);
    S3_m(m-5:m) = fliplr(-S3_U);


    H = h*H;
    HI = inv(H);

    % Fourth derivative, 1th order accurate at first 8 boundary points (still
    % yield 5th order convergence if stable: for example u_tt = -u_xxxx

    m4 = 7/240;
    m3 = -2/5;
    m2 = 169/60;
    m1 = -122/15;
    m0 = 91/8;

    M4 = m4*(diag(ones(m-4,1),4)+diag(ones(m-4,1),-4))+m3*(diag(ones(m-3,1),3)+diag(ones(m-3,1),-3))+m2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2))+m1*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1))+m0*diag(ones(m,1),0);

    %M4 = (-1/6*(diag(ones(m-3,1),3)+diag(ones(m-3,1),-3) ) + 2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2)) -13/2*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1)) + 28/3*diag(ones(m,1),0));

    M4_U = [
        0.1009e4/0.192e3 -0.7657e4/0.480e3 0.9307e4/0.480e3 -0.509e3/0.40e2 0.4621e4/0.960e3 -0.25e2/0.32e2;
        -0.7657e4/0.480e3 0.49513e5/0.960e3 -0.4007e4/0.60e2 0.21799e5/0.480e3 -0.8171e4/0.480e3 0.2657e4/0.960e3;
        0.9307e4/0.480e3 -0.4007e4/0.60e2 0.1399e4/0.15e2 -0.2721e4/0.40e2 0.12703e5/0.480e3 -0.521e3/0.120e3;
        -0.509e3/0.40e2 0.21799e5/0.480e3 -0.2721e4/0.40e2 0.3349e4/0.60e2 -0.389e3/0.15e2 0.559e3/0.96e2;
        0.4621e4/0.960e3 -0.8171e4/0.480e3 0.12703e5/0.480e3 -0.389e3/0.15e2 0.17857e5/0.960e3 -0.1499e4/0.160e3;
        -0.25e2/0.32e2 0.2657e4/0.960e3 -0.521e3/0.120e3 0.559e3/0.96e2 -0.1499e4/0.160e3 0.2225e4/0.192e3;
    ];

    M4(1:6,1:6) = M4_U;

    M4(m-5:m,m-5:m) = flipud( fliplr( M4_U ) );
    M4 = M4/h^3;

    D4 = HI*(M4-e_1*S3_1+e_m*S3_m  + S_1'*S2_1-S_m'*S2_m);
end
