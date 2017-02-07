function [H, HI, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = d4_variable_4_min_boundary_points(m,h)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% 4:de ordn. SBP Finita differens         %%%
    %%% operatorer framtagna av Mark Carpenter  %%%
    %%%                                         %%%
    %%% H           (Normen)                    %%%
    %%% D1=H^(-1)Q  (approx f?rsta derivatan)   %%%
    %%% D2          (approx andra derivatan)    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %H?r med endast 4 randpunkter


    BP = 4;
    if(m<2*BP)
        error(['Operator requires at least ' num2str(2*BP) ' grid points']);
    end


    % Norm
    Hv = ones(m,1);
    Hv(1:4) = [17/48 59/48 43/48 49/48];
    Hv(m-3:m) = rot90(Hv(1:4),2);
    Hv = h*Hv;
    H = spdiag(Hv, 0);
    HI = spdiag(1./Hv, 0);


    % Boundary operators
    e_l = sparse(m,1);
    e_l(1) = 1;
    e_r = rot90(e_l, 2);

    d1_l = sparse(m,1);
    d1_l(1:4) = 1/h*[-11/6 3 -3/2 1/3];
    d1_r = -rot90(d1_l, 2);

    d2_l = sparse(m,1);
    d2_l(1:4) = 1/h^2*[2 -5 4 -1];
    d2_r = rot90(d2_l, 2);

    d3_l = sparse(m,1);
    d3_l(1:4) = 1/h^3*[-1 3 -3 1];
    d3_r = -rot90(d3_l, 2);


    % First derivative
    stencil = [1/12 -2/3 0 2/3 -1/12];
    diags = [-1 0 1];

    Q_U = [
        0 0.59e2/0.96e2 -0.1e1/0.12e2 -0.1e1/0.32e2;
         -0.59e2/0.96e2 0 0.59e2/0.96e2 0;
         0.1e1/0.12e2 -0.59e2/0.96e2 0 0.59e2/0.96e2;
         0.1e1/0.32e2 0 -0.59e2/0.96e2 0;
    ];

    Q = stripeMatrix(stencil, diags, m);
    Q(1:4,1:4)=Q_U;
    Q(m-3:m,m-3:m) = -rot90(Q_U, 2);

    D1 = HI*(Q - 1/2*e_l*e_l' + 1/2*e_r*e_r');

    % Fourth derivative
    stencil = [-1/6, 2, -13/2, 28/3, -13/2, 2, -1/6];
    diags = -3:3;
    M4 = stripeMatrix(stencil, diags, m);

    M4_U=[
        0.8e1/0.3e1 -0.37e2/0.6e1 0.13e2/0.3e1 -0.5e1/0.6e1;
        -0.37e2/0.6e1 0.47e2/0.3e1 -13 0.11e2/0.3e1;
        0.13e2/0.3e1 -13 0.44e2/0.3e1 -0.47e2/0.6e1;
        -0.5e1/0.6e1 0.11e2/0.3e1 -0.47e2/0.6e1 0.29e2/0.3e1;
    ];


    M4(1:4,1:4) = M4_U;
    M4(m-3:m,m-3:m) = rot90(M4_U, 2);
    M4 = 1/h^3*M4;

    D4=HI*(M4 - e_l*d3_l'+e_r*d3_r' + d1_l*d2_l'-d1_r*d2_r');
end