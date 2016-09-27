function [H, HI, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = d4_variable_6_min_boundary_points(m,h)
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

    BP = 6;
    if(m<2*BP)
        error(['Operator requires at least ' num2str(2*BP) ' grid points']);
    end

    % Norm
    Hv = ones(m,1);
    Hv(1:6) = [13649/43200,12013/8640,2711/4320,5359/4320,7877/8640, 43801/43200];
    Hv(m-5:m) = rot90(Hv(1:6),2);
    Hv = h*Hv;
    H = spdiag(Hv, 0);
    HI = spdiag(1./Hv, 0);


    % Boundary operators
    e_l = sparse(m,1);
    e_l(1) = 1;
    e_r = rot90(e_l, 2);

    d1_l = sparse(m,1);
    d1_l(1:5) = [-25/12, 4, -3, 4/3, -1/4]/h;
    d1_r = -rot90(d1_l, 2);

    d2_l = sparse(m,1);
    d2_l(1:5) = [0.35e2/0.12e2 -0.26e2/0.3e1 0.19e2/0.2e1 -0.14e2/0.3e1 0.11e2/0.12e2;]/h^2;
    d2_r = rot90(d2_l, 2);

    d3_l = sparse(m,1);
    d3_l(1:5) = [-0.5e1/0.2e1 9 -12 7 -0.3e1/0.2e1;]/h^3;
    d3_r = -rot90(d3_l, 2);


    % Fourth derivative, 1th order accurate at first 8 boundary points (still
    % yield 5th order convergence if stable: for example u_tt=-u_xxxx

    stencil = [7/240, -2/5, 169/60, -122/15, 91/8, -122/15, 169/60, -2/5, 7/240];
    diags = -4:4;
    M4 = stripeMatrix(stencil, diags, m);

    M4_U=[
        0.3504379e7/0.907200e6 -0.4613983e7/0.453600e6 0.4260437e7/0.453600e6 -0.418577e6/0.113400e6 0.524579e6/0.907200e6 0.535e3/0.18144e5;
        -0.4613983e7/0.453600e6 0.5186159e7/0.181440e6 -0.81121e5/0.2835e4 0.218845e6/0.18144e5 -0.159169e6/0.90720e5 -0.94669e5/0.907200e6;
        0.4260437e7/0.453600e6 -0.81121e5/0.2835e4 0.147695e6/0.4536e4 -0.384457e6/0.22680e5 0.339653e6/0.90720e5 -0.18233e5/0.113400e6;
        -0.418577e6/0.113400e6 0.218845e6/0.18144e5 -0.384457e6/0.22680e5 0.65207e5/0.4536e4 -0.22762e5/0.2835e4 0.1181753e7/0.453600e6;
        0.524579e6/0.907200e6 -0.159169e6/0.90720e5 0.339653e6/0.90720e5 -0.22762e5/0.2835e4 0.2006171e7/0.181440e6 -0.3647647e7/0.453600e6;
        0.535e3/0.18144e5 -0.94669e5/0.907200e6 -0.18233e5/0.113400e6 0.1181753e7/0.453600e6 -0.3647647e7/0.453600e6 0.10305271e8/0.907200e6;
    ];

    M4(1:6,1:6) = M4_U;
    M4(m-5:m,m-5:m) = rot90(M4_U, 2);
    M4 = 1/h^3*M4;

    D4=HI*(M4 - e_l*d3_l'+e_r*d3_r' + d1_l*d2_l'-d1_r*d2_r');
end
