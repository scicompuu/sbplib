% Returns D2 as a function handle
function [H, HI, D1, D2, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = d4_variable_2(m,h)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% 4:de ordn. SBP Finita differens         %%%
    %%% operatorer framtagna av Ken Mattsson    %%%
    %%%                                         %%%
    %%% 6 randpunkter, diagonal norm            %%%
    %%%                                         %%%
    %%% Datum: 2013-11-11                       %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    BP = 2;
    if(m < 2*BP)
        error('Operator requires at least %d grid points', 2*BP);
    end

    % Norm
    Hv = ones(m,1);
    Hv(1) = 1/2;
    Hv(m) = 1/2;
    Hv = h*Hv;
    H = spdiag(Hv, 0);
    HI = spdiag(1./Hv, 0);

    % Boundary operators
    e_l = sparse(m,1);
    e_l(1) = 1;
    e_r = rot90(e_l, 2);

    d1_l = sparse(m,1);
    d1_l(1:3) = 1/h*[-3/2 2 -1/2];
    d1_r = -rot90(d1_l, 2);

    d2_l = sparse(m,1);
    d2_l(1:3) = 1/h^2*[1 -2 1];
    d2_r = rot90(d2_l, 2);

    d3_l = sparse(m,1);
    d3_l(1:4) = 1/h^3*[-1 3 -3 1];
    d3_r = -rot90(d3_l, 2);


    % First derivative SBP operator, 1st order accurate at first 6 boundary points
    stencil = [-1/2, 0, 1/2];
    diags = [-1 0 1];
    Q = stripeMatrix(stencil, diags, m);

    D1 = HI*(Q - 1/2*e_l*e_l' + 1/2*e_r*e_r');

    % Second derivative, 1st order accurate at first boundary points
    M = sparse(m,m);

    scheme_width = 3;
    scheme_radius = (scheme_width-1)/2;
    r = (1+scheme_radius):(m-scheme_radius);

    function D2 = D2_fun(c)
        Mm1 = -c(r-1)/2 - c(r)/2;
        M0  =  c(r-1)/2 + c(r)   + c(r+1)/2;
        Mp1 =            -c(r)/2 - c(r+1)/2;

        M(r,:) = spdiags([Mm1 M0 Mp1],0:2*scheme_radius,length(r),m);

        M(1:2,1:2) = [c(1)/2 + c(2)/2 -c(1)/2 - c(2)/2; -c(1)/2 - c(2)/2 c(1)/2 + c(2) + c(3)/2;];
        M(m-1:m,m-1:m) = [c(m-2)/2 + c(m-1) + c(m)/2 -c(m-1)/2 - c(m)/2; -c(m-1)/2 - c(m)/2 c(m-1)/2 + c(m)/2;];
        M = 1/h*M;

        D2 = HI*(-M - c(1)*e_l*d1_l' + c(m)*e_r*d1_r');
    end
    D2 = @D2_fun;

    % Fourth derivative, 0th order accurate at first 6 boundary points
    stencil = [1, -4, 6, -4, 1];
    diags = -2:2;
    M4 = stripeMatrix(stencil, diags, m);

    M4_U = [
         0.13e2/0.10e2 -0.12e2/0.5e1   0.9e1/0.10e2   0.1e1/0.5e1;
        -0.12e2/0.5e1   0.26e2/0.5e1  -0.16e2/0.5e1   0.2e1/0.5e1;
         0.9e1/0.10e2  -0.16e2/0.5e1   0.47e2/0.10e2 -0.17e2/0.5e1;
         0.1e1/0.5e1    0.2e1/0.5e1   -0.17e2/0.5e1   0.29e2/0.5e1;
    ];

    M4(1:4,1:4) = M4_U;
    M4(m-3:m,m-3:m) = rot90(M4_U, 2);
    M4 = 1/h^3*M4;

    D4=HI*(M4 - e_l*d3_l'+e_r*d3_r' + d1_l*d2_l'-d1_r*d2_r');
end
