function [H, HI, Dp, Dm, e_1, e_m] = d1_upwind_4(m,h)

    if(m<8)
        error('Operator requires at least 8 grid points');
    end

    Hv=ones(m,1);
    Hv(1:4)=[49/144 61/48 41/48 149/144];
    Hv(m-3:m)=rot90(Hv(1:4),2);
    Hv = Hv*h;
    H = spdiag(Hv,0);
    HI = spdiag(1./Hv,0);

    q_diags   = [-1, 0, 1, 2 3];
    q_stencil = [-1/4 -5/6 3/2 -1/2 1/12];
    Qp = stripeMatrix(q_stencil, q_diags,m);

    Q_U = [
        -0.1e1/0.48e2    0.205e3/0.288e3 -0.29e2/0.144e3 0.1e1/0.96e2;
        -0.169e3/0.288e3 -0.11e2/0.48e2  0.33e2/0.32e2   -0.43e2/0.144e3;
        0.11e2/0.144e3   -0.13e2/0.32e2  -0.29e2/0.48e2  0.389e3/0.288e3;
        0.1e1/0.32e2     -0.11e2/0.144e3 -0.65e2/0.288e3 -0.13e2/0.16e2;
    ];

    Qp(1:4,1:4)=Q_U;
    Qp(m-3:m,m-3:m)=rot90(Q_U, 2)'; %%% This is different from standard SBP

    Qm=-Qp';

    e_1=sparse(m,1);e_1(1)=1;
    e_m=sparse(m,1);e_m(m)=1;

    Dp=HI*(Qp-1/2*(e_1*e_1')+1/2*(e_m*e_m')) ;

    Dm=HI*(Qm-1/2*(e_1*e_1')+1/2*(e_m*e_m')) ;
end
