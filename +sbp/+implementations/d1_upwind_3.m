function [H, HI, Dp, Dm, e_1, e_m] = d1_upwind_3(m,h)
    
    if(m<6)
        error('Operator requires at least 6 grid points');
    end

    Hv = ones(m,1);
    Hv(1:3) = [3/8; 7/6; 23/24];
    Hv(m-2:m) = rot90(Hv(1:3),2);
    Hv = Hv*h;
    H = spdiag(Hv,0);
    HI = spdiag(1./Hv,0);

    q_diags   = [-1, 0, 1, 2];
    q_stencil = [-1/3 -1/2 1 -1/6];
    Qp = stripeMatrix(q_stencil, q_diags,m);

    Q_U = [
         -1/24  17/24   -1/6;
        -13/24   -1/4  23/24;
          1/12 -11/24 -11/24;
    ];
    Qp(1:3,1:3)=Q_U;
    Qp(m-2:m,m-2:m)=rot90(Q_U,2)'; %%% This is different from standard SBP

    Qm=-Qp';

    e_1=sparse(m,1);
    e_1(1)=1;
    e_m=sparse(m,1);
    e_m(m)=1;

    Dp=HI*(Qp-1/2*(e_1*e_1')+1/2*(e_m*e_m')) ;

    Dm=HI*(Qm-1/2*(e_1*e_1')+1/2*(e_m*e_m')) ;
end

