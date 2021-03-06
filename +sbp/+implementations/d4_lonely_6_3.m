function [H, HI, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = d4_variable_6_3(m,h)
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

    % H?r med 7 RP ist?llet f?r 8 f?r D4 operatorn, dock samma randderivator
    % Denna ?r noggrannare, och har 2a ordningens randdslutning och b?r ge 6te
    % ordningens konvergens. Hade 2 fria parametrar att optimera

    % Norm
    Hv = ones(m,1);
    Hv(1:7) = [0.414837907e9/0.1191965760e10, 0.475278367e9/0.397321920e9, 0.13872751e8/0.12416310e8, 0.346739027e9/0.595982880e9, 0.560227469e9/0.397321920e9, 0.322971631e9/0.397321920e9, 0.616122491e9/0.595982880e9];
    Hv(m-6:m) = rot90(Hv(1:7),2);
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


    % Fourth derivative, 1th order accurate at first 8 boundary points
    stencil = [7/240, -2/5, 169/60, -122/15, 91/8, -122/15, 169/60, -2/5, 7/240];
    diags = -4:4;
    M4 = stripeMatrix(stencil, diags, m);

    M4_U = [
        0.1399708478939e13/0.263487168000e12 -0.13482796013041e14/0.834376032000e12 0.344344095859e12/0.17565811200e11 -0.3166261424681e13/0.250312809600e12 0.1508605165681e13/0.333750412800e12 -0.486270829441e12/0.834376032000e12 -0.221976356359e12/0.5006256192000e13;
        -0.13482796013041e14/0.834376032000e12 0.7260475818391e13/0.139062672000e12 -0.27224036353e11/0.406022400e9 0.1847477458951e13/0.41718801600e11 -0.848984558161e12/0.55625068800e11 0.247494925991e12/0.139062672000e12 0.165585445559e12/0.834376032000e12;
        0.344344095859e12/0.17565811200e11 -0.27224036353e11/0.406022400e9 0.2044938640393e13/0.22250027520e11 -0.1071086785417e13/0.16687520640e11 0.502199537033e12/0.22250027520e11 -0.143589154441e12/0.55625068800e11 -0.88181965559e11/0.333750412800e12;
        -0.3166261424681e13/0.250312809600e12 0.1847477458951e13/0.41718801600e11 -0.1071086785417e13/0.16687520640e11 0.628860435593e12/0.12515640480e11 -0.73736245829e11/0.3337504128e10 0.195760572271e12/0.41718801600e11 -0.81156046361e11/0.250312809600e12;
        0.1508605165681e13/0.333750412800e12 -0.848984558161e12/0.55625068800e11 0.502199537033e12/0.22250027520e11 -0.73736245829e11/0.3337504128e10 0.76725285869e11/0.4450005504e10 -0.3912429433e10/0.406022400e9 0.53227370659e11/0.17565811200e11;
        -0.486270829441e12/0.834376032000e12 0.247494925991e12/0.139062672000e12 -0.143589154441e12/0.55625068800e11 0.195760572271e12/0.41718801600e11 -0.3912429433e10/0.406022400e9 0.1699707221791e13/0.139062672000e12 -0.6959018412841e13/0.834376032000e12;
        -0.221976356359e12/0.5006256192000e13 0.165585445559e12/0.834376032000e12 -0.88181965559e11/0.333750412800e12 -0.81156046361e11/0.250312809600e12 0.53227370659e11/0.17565811200e11 -0.6959018412841e13/0.834376032000e12 0.3012195053939e13/0.263487168000e12;
    ];

    M4(1:7,1:7) = M4_U;
    M4(m-6:m,m-6:m) = rot90(M4_U, 2);
    M4 = 1/h^3*M4;

    D4=HI*(M4 - e_l*d3_l'+e_r*d3_r' + d1_l*d2_l'-d1_r*d2_r');
end
