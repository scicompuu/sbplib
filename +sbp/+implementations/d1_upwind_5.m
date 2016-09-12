function [H, HI, Dp, Dm, e_1, e_m] = d1_upwind_5(m,h)
    
    if(m<8)
        error('Operator requires at least 8 grid points');
    end

    Hv=ones(m,1);
    Hv(1:4)=[0.251e3/0.720e3; 0.299e3/0.240e3; 0.211e3/0.240e3; 0.739e3/0.720e3;];
    Hv(m-3:m)=rot90(Hv(1:4),2);
    Hv = Hv*h;
    H = spdiag(Hv,0);
    HI = spdiag(1./Hv,0);

    q_diags   = [-2 -1 0 1 2 3];
    q_stencil = [1/20 -1/2 -1/3 +1 -1/4 +1/30];
    Qp = stripeMatrix(q_stencil, q_diags,m);

    Q_U = [
        -0.1e1/0.120e3 0.941e3/0.1440e4 -0.47e2/0.360e3 -0.7e1/0.480e3;
        -0.869e3/0.1440e4 -0.11e2/0.120e3 0.25e2/0.32e2 -0.43e2/0.360e3;
        0.29e2/0.360e3 -0.17e2/0.32e2 -0.29e2/0.120e3 0.1309e4/0.1440e4;
        0.1e1/0.32e2 -0.11e2/0.360e3 -0.661e3/0.1440e4 -0.13e2/0.40e2;
    ];

    Qp(1:4,1:4)=Q_U;
    Qp(m-3:m,m-3:m)=rot90( Q_U(1:4,1:4) ,2 )'; %%% This is different from standard SBP

    Qm=-Qp';

    e_1=sparse(m,1);e_1(1)=1;
    e_m=sparse(m,1);e_m(m)=1;

    Dp=HI*(Qp-1/2*(e_1*e_1')+1/2*(e_m*e_m')) ;

    Dm=HI*(Qm-1/2*(e_1*e_1')+1/2*(e_m*e_m')) ;
end
