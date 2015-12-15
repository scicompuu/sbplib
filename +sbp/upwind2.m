function [H, HI, Dp, Dm, e_1, e_m] = upwind2(m,h)
    Hv=ones(m,1);
    Hv(1:2)=[0.1e1/0.4e1; 0.5e1/0.4e1;];
    Hv(m-1:m)=rot90(Hv(1:2),2);
    Hv = Hv*h;
    H = spdiag(Hv,0);
    HI = spdiag(1./Hv,0);

    q_diags   = [0 1 2];
    q_stencil = [-3/2 +2 -1/2];
    Qp = stripeMatrix(q_stencil, q_diags,m);

    Q_U = [
        -0.1e1/0.4e1 0.5e1/0.4e1;
        -0.1e1/0.4e1 -0.5e1/0.4e1;
    ];

    Qp(1:2,1:2)=Q_U;
    Qp(m-1:m,m-1:m)=rot90(Q_U,2)'; %%% This is different from standard SBP

    Qm=-Qp';

    e_1=zeros(m,1);e_1(1)=1;
    e_m=zeros(m,1);e_m(m)=1;

    Dp=HI*(Qp-1/2*e_1*e_1'+1/2*e_m*e_m') ;

    Dm=HI*(Qm-1/2*e_1*e_1'+1/2*e_m*e_m') ;
end



% KOLLA NAMN PÅ ALLA FUNKTIONER
% TESTA ALLA FUNKTIONER