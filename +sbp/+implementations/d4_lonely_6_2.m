function [H, HI, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = d4_variable_6_2(m,h)
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

    BP = 6;
    if(m<2*BP)
        error(['Operator requires at least ' num2str(2*BP) ' grid points']);
    end

    % Norm
    Hv = ones(m,1);
    Hv(1:6) = [0.181e3/0.576e3, 0.1343e4/0.960e3, 0.293e3/0.480e3, 0.1811e4/0.1440e4, 0.289e3/0.320e3, 0.65e2/0.64e2];
    Hv(m-5:m) = rot90(Hv(1:6),2);
    Hv = h*Hv;
    H = spdiag(Hv, 0);
    HI = spdiag(1./Hv, 0);


    % Boundary operators
    e_l = sparse(m,1);
    e_l(1) = 1;
    e_r = rot90(e_l, 2);

    d1_l = sparse(m,1);
    d1_l(1:6) = [-0.137e3/0.60e2 5 -5 0.10e2/0.3e1 -0.5e1/0.4e1 0.1e1/0.5e1;]/h;
    d1_r = -rot90(d1_l, 2);

    d2_l = sparse(m,1);
    d2_l(1:6) = [0.15e2/0.4e1 -0.77e2/0.6e1 0.107e3/0.6e1 -13 0.61e2/0.12e2 -0.5e1/0.6e1;]/h^2;
    d2_r = rot90(d2_l, 2);

    d3_l = sparse(m,1);
    d3_l(1:6) = [-0.17e2/0.4e1 0.71e2/0.4e1 -0.59e2/0.2e1 0.49e2/0.2e1 -0.41e2/0.4e1 0.7e1/0.4e1;]/h^3;
    d3_r = -rot90(d3_l, 2);


    % Fourth derivative, 1th order accurate at first 8 boundary points (still
    % yield 5th order convergence if stable: for example u_tt = -u_xxxx
    stencil = [7/240, -2/5, 169/60, -122/15, 91/8, -122/15, 169/60, -2/5, 7/240];
    diags = -4:4;
    M4 = stripeMatrix(stencil, diags, m);

    M4_U = [
        0.1009e4/0.192e3 -0.7657e4/0.480e3 0.9307e4/0.480e3 -0.509e3/0.40e2 0.4621e4/0.960e3 -0.25e2/0.32e2;
        -0.7657e4/0.480e3 0.49513e5/0.960e3 -0.4007e4/0.60e2 0.21799e5/0.480e3 -0.8171e4/0.480e3 0.2657e4/0.960e3;
        0.9307e4/0.480e3 -0.4007e4/0.60e2 0.1399e4/0.15e2 -0.2721e4/0.40e2 0.12703e5/0.480e3 -0.521e3/0.120e3;
        -0.509e3/0.40e2 0.21799e5/0.480e3 -0.2721e4/0.40e2 0.3349e4/0.60e2 -0.389e3/0.15e2 0.559e3/0.96e2;
        0.4621e4/0.960e3 -0.8171e4/0.480e3 0.12703e5/0.480e3 -0.389e3/0.15e2 0.17857e5/0.960e3 -0.1499e4/0.160e3;
        -0.25e2/0.32e2 0.2657e4/0.960e3 -0.521e3/0.120e3 0.559e3/0.96e2 -0.1499e4/0.160e3 0.2225e4/0.192e3;
    ];


    M4(1:6,1:6) = M4_U;
    M4(m-5:m,m-5:m) = rot90(M4_U, 2);
    M4 = 1/h^3*M4;

    D4=HI*(M4 - e_l*d3_l'+e_r*d3_r' + d1_l*d2_l'-d1_r*d2_r');
end
