function [H, HI, Dp, Dm, e_1, e_m] = d1_upwind_8(m,h)
    
    BP = 8;
    if(m<2*BP)
        error(['Operator requires at least ' num2str(2*BP) ' grid points']);
    end

    Hv=ones(m,1);
    Hv(1:8) = [0.7489399e7/0.25401600e8; 0.5537831e7/0.3628800e7; 0.103373e6/0.403200e6; 0.261259e6/0.145152e6; 0.298231e6/0.725760e6; 0.515917e6/0.403200e6; 0.3349159e7/0.3628800e7; 0.25639991e8/0.25401600e8];
    Hv(m-7:m)=rot90(Hv(1:8),2);
    Hv=Hv*h;

    H = spdiag(Hv,0);
    HI = spdiag(1./Hv,0);

    q_diags   = [-3 -2 -1 0 1 2 3 4 5];
    q_stencil = [-1/168 +1/14 -1/2 -9/20 +5/4 -1/2 +1/6 -1/28 +1/280];
    Qp = stripeMatrix(q_stencil, q_diags,m);


    Q_U =[
        -0.16683e5/0.63018560e8 0.29345822969987e14/0.43354248537600e14 -0.2734625426591e13/0.40644608004000e14 -0.113480208109603e15/0.780376473676800e15 -0.830250230261e12/0.26012549122560e14 0.2500519492033e13/0.32515686403200e14 0.2235718279643e13/0.390188236838400e15 -0.388481888477e12/0.26543417472000e14;
        -0.29227665839987e14/0.43354248537600e14 -0.493793e6/0.63018560e8 0.8302717120817e13/0.26543417472000e14 0.3739408501537e13/0.9290196115200e13 0.2684481534461e13/0.13935294172800e14 -0.4450185662513e13/0.18580392230400e14 -0.1221838279381e13/0.37160784460800e14 0.90595000956023e14/0.1950941184192000e16;
        0.2505689537591e13/0.40644608004000e14 -0.7312922392817e13/0.26543417472000e14 -0.69332623e8/0.1323389760e10 0.10994933811709e14/0.18580392230400e14 -0.9270952411151e13/0.18580392230400e14 0.3191238635141e13/0.20644880256000e14 0.4442211176987e13/0.92901961152000e14 -0.940661365031e12/0.32515686403200e14;
        0.118016946570403e15/0.780376473676800e15 -0.4173878828737e13/0.9290196115200e13 -0.7990503962509e13/0.18580392230400e14 -0.207799621e9/0.1323389760e10 0.2044021560341e13/0.2477385630720e13 0.511197701761e12/0.18580392230400e14 0.1237681717213e13/0.13935294172800e14 -0.7784834666617e13/0.130062745612800e15;
        0.68609076271e11/0.2364777192960e13 -0.2235651762161e13/0.13935294172800e14 0.6527681584751e13/0.18580392230400e14 -0.1115980068821e13/0.2477385630720e13 -0.55386253e8/0.189055680e9 0.3208334350649e13/0.3716078446080e13 -0.407569013461e12/0.844563283200e12 0.136474842626653e15/0.780376473676800e15;
        -0.2487637785013e13/0.32515686403200e14 0.4244231077313e13/0.18580392230400e14 -0.1550378843141e13/0.20644880256000e14 -0.5726967564961e13/0.18580392230400e14 -0.1017898941929e13/0.3716078446080e13 -0.526653889e9/0.1323389760e10 0.45241297077547e14/0.37160784460800e14 -0.2279608411897e13/0.5080576000500e13;
        -0.2164019088443e13/0.390188236838400e15 0.1263196075861e13/0.37160784460800e14 -0.6600697610987e13/0.92901961152000e14 0.556610591687e12/0.13935294172800e14 0.926842346471e12/0.9290196115200e13 -0.18757693936747e14/0.37160784460800e14 -0.584765899e9/0.1323389760e10 0.204646287449e12/0.168431424000e12;
        0.387091928477e12/0.26543417472000e14 -0.90231551688023e14/0.1950941184192000e16 0.1032404418251e13/0.32515686403200e14 0.3502353445417e13/0.130062745612800e15 -0.15385068876253e14/0.780376473676800e15 0.262499068919e12/0.10161152001000e14 -0.867004691939e12/0.1852745664000e13 -0.85017967e8/0.189055680e9;
    ];

    Qp(1:8,1:8)=Q_U;
    Qp(m-7:m,m-7:m)=rot90(Q_U,2)'; %%% This is different from standard SBP

    Qm=-Qp';

    e_1=sparse(m,1);e_1(1)=1;
    e_m=sparse(m,1);e_m(m)=1;

    Dp=HI*(Qp-1/2*e_1*e_1'+1/2*e_m*e_m') ;

    Dm=HI*(Qm-1/2*e_1*e_1'+1/2*e_m*e_m') ;
end
