function [H, HI, D1, D2, e_l, e_r, d1_l, d1_r] = d2_variable_4(m,h)

    BP = 6;
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




    S = d1_l*d1_l' + d1_r*d1_r';

    stencil = [1/12 -2/3 0 2/3 -1/12];
    diags = -2:2;

    Q_U = [
        0 0.59e2/0.96e2 -0.1e1/0.12e2 -0.1e1/0.32e2;
        -0.59e2/0.96e2 0 0.59e2/0.96e2 0;
        0.1e1/0.12e2 -0.59e2/0.96e2 0 0.59e2/0.96e2;
        0.1e1/0.32e2 0 -0.59e2/0.96e2 0;
    ];

    Q = stripeMatrix(stencil, diags, m);
    Q(1:4,1:4) = Q_U;
    Q(m-3:m,m-3:m) = -rot90(Q_U, 2);

    D1 = HI*(Q - 1/2*e_l*e_l' + 1/2*e_r*e_r');


    N = m;
    function D2 = D2_fun(c)
        M = 78+(N-12)*5;
        %h = 1/(N-1);


        U = [
            48/(17)*(0.12e2/0.17e2 * c(1) + 0.59e2/0.192e3 * c(2) + 0.27010400129e11/0.345067064608e12 * c(3) + 0.69462376031e11/0.2070402387648e13 * c(4)) 48/(17)*(-0.59e2/0.68e2 * c(1) - 0.6025413881e10/0.21126554976e11 * c(3) - 0.537416663e9/0.7042184992e10 * c(4)) 48/(17)*(0.2e1/0.17e2 * c(1) - 0.59e2/0.192e3 * c(2) + 0.213318005e9/0.16049630912e11 * c(4) + 0.2083938599e10/0.8024815456e10 * c(3)) 48/(17)*(0.3e1/0.68e2 * c(1) - 0.1244724001e10/0.21126554976e11 * c(3) + 0.752806667e9/0.21126554976e11 * c(4)) 48/(17)*(0.49579087e8/0.10149031312e11 * c(3) - 0.49579087e8/0.10149031312e11 * c(4)) 48/(17)*(-c(4)/0.784e3 + c(3)/0.784e3);...
            48/(59)*(-0.59e2/0.68e2 * c(1) - 0.6025413881e10/0.21126554976e11 * c(3) - 0.537416663e9/0.7042184992e10 * c(4)) 48/(59)*(0.3481e4/0.3264e4 * c(1) + 0.9258282831623875e16/0.7669235228057664e16 * c(3) + 0.236024329996203e15/0.1278205871342944e16 * c(4)) 48/(59)*(-0.59e2/0.408e3 * c(1) - 0.29294615794607e14/0.29725717938208e14 * c(3) - 0.2944673881023e13/0.29725717938208e14 * c(4)) 48/(59)*(-0.59e2/0.1088e4 * c(1) + 0.260297319232891e15/0.2556411742685888e16 * c(3) - 0.60834186813841e14/0.1278205871342944e16 * c(4)) 48/(59)*(-0.1328188692663e13/0.37594290333616e14 * c(3) + 0.1328188692663e13/0.37594290333616e14 * c(4)) 48/(59)*(-0.8673e4/0.2904112e7 * c(3) + 0.8673e4/0.2904112e7 * c(4));...
            48/(43)*(0.2e1/0.17e2 * c(1) - 0.59e2/0.192e3 * c(2) + 0.213318005e9/0.16049630912e11 * c(4) + 0.2083938599e10/0.8024815456e10 * c(3)) 48/(43)*(-0.59e2/0.408e3 * c(1) - 0.29294615794607e14/0.29725717938208e14 * c(3) - 0.2944673881023e13/0.29725717938208e14 * c(4)) 48/(43)*(c(1)/0.51e2 + 0.59e2/0.192e3 * c(2) + 0.13777050223300597e17/0.26218083221499456e17 * c(4) + 0.564461e6/0.13384296e8 * c(5) + 0.378288882302546512209e21/0.270764341349677687456e21 * c(3)) 48/(43)*(c(1)/0.136e3 - 0.125059e6/0.743572e6 * c(5) - 0.4836340090442187227e19/0.5525802884687299744e19 * c(3) - 0.17220493277981e14/0.89177153814624e14 * c(4)) 48/(43)*(-0.10532412077335e14/0.42840005263888e14 * c(4) + 0.1613976761032884305e19/0.7963657098519931984e19 * c(3) + 0.564461e6/0.4461432e7 * c(5)) 48/(43)*(-0.960119e6/0.1280713392e10 * c(4) - 0.3391e4/0.6692148e7 * c(5) + 0.33235054191e11/0.26452850508784e14 * c(3));...
            48/(49)*(0.3e1/0.68e2 * c(1) - 0.1244724001e10/0.21126554976e11 * c(3) + 0.752806667e9/0.21126554976e11 * c(4)) 48/(49)*(-0.59e2/0.1088e4 * c(1) + 0.260297319232891e15/0.2556411742685888e16 * c(3) - 0.60834186813841e14/0.1278205871342944e16 * c(4)) 48/(49)*(c(1)/0.136e3 - 0.125059e6/0.743572e6 * c(5) - 0.4836340090442187227e19/0.5525802884687299744e19 * c(3) - 0.17220493277981e14/0.89177153814624e14 * c(4)) 48/(49)*(0.3e1/0.1088e4 * c(1) + 0.507284006600757858213e21/0.475219048083107777984e21 * c(3) + 0.1869103e7/0.2230716e7 * c(5) + c(6)/0.24e2 + 0.1950062198436997e16/0.3834617614028832e16 * c(4)) 48/(49)*(-0.4959271814984644613e19/0.20965546238960637264e20 * c(3) - c(6)/0.6e1 - 0.15998714909649e14/0.37594290333616e14 * c(4) - 0.375177e6/0.743572e6 * c(5)) 48/(49)*(-0.368395e6/0.2230716e7 * c(5) + 0.752806667e9/0.539854092016e12 * c(3) + 0.1063649e7/0.8712336e7 * c(4) + c(6)/0.8e1);...
            0.49579087e8/0.10149031312e11 * c(3) - 0.49579087e8/0.10149031312e11 * c(4) -0.1328188692663e13/0.37594290333616e14 * c(3) + 0.1328188692663e13/0.37594290333616e14 * c(4) -0.10532412077335e14/0.42840005263888e14 * c(4) + 0.1613976761032884305e19/0.7963657098519931984e19 * c(3) + 0.564461e6/0.4461432e7 * c(5) -0.4959271814984644613e19/0.20965546238960637264e20 * c(3) - c(6)/0.6e1 - 0.15998714909649e14/0.37594290333616e14 * c(4) - 0.375177e6/0.743572e6 * c(5) 0.8386761355510099813e19/0.128413970713633903242e21 * c(3) + 0.2224717261773437e16/0.2763180339520776e16 * c(4) + 0.5e1/0.6e1 * c(6) + c(7)/0.24e2 + 0.280535e6/0.371786e6 * c(5) -0.35039615e8/0.213452232e9 * c(4) - c(7)/0.6e1 - 0.13091810925e11/0.13226425254392e14 * c(3) - 0.1118749e7/0.2230716e7 * c(5) - c(6)/0.2e1;...
            -c(4)/0.784e3 + c(3)/0.784e3 -0.8673e4/0.2904112e7 * c(3) + 0.8673e4/0.2904112e7 * c(4) -0.960119e6/0.1280713392e10 * c(4) - 0.3391e4/0.6692148e7 * c(5) + 0.33235054191e11/0.26452850508784e14 * c(3) -0.368395e6/0.2230716e7 * c(5) + 0.752806667e9/0.539854092016e12 * c(3) + 0.1063649e7/0.8712336e7 * c(4) + c(6)/0.8e1 -0.35039615e8/0.213452232e9 * c(4) - c(7)/0.6e1 - 0.13091810925e11/0.13226425254392e14 * c(3) - 0.1118749e7/0.2230716e7 * c(5) - c(6)/0.2e1 0.3290636e7/0.80044587e8 * c(4) + 0.5580181e7/0.6692148e7 * c(5) + 0.5e1/0.6e1 * c(7) + c(8)/0.24e2 + 0.660204843e9/0.13226425254392e14 * c(3) + 0.3e1/0.4e1 * c(6)
        ];


        L = [
            c(N-7)/0.24e2 + 0.5e1/0.6e1 * c(N-6) + 0.5580181e7/0.6692148e7 * c(N-4) + 0.4887707739997e13/0.119037827289528e15 * c(N-3) + 0.3e1/0.4e1 * c(N-5) + 0.660204843e9/0.13226425254392e14 * c(N-2) + 0.660204843e9/0.13226425254392e14 * c(N-1) -c(N-6)/0.6e1 - 0.1618585929605e13/0.9919818940794e13 * c(N-3) - c(N-5)/0.2e1 - 0.1118749e7/0.2230716e7 * c(N-4) - 0.13091810925e11/0.13226425254392e14 * c(N-2) - 0.13091810925e11/0.13226425254392e14 * c(N-1) -0.368395e6/0.2230716e7 * c(N-4) + c(N-5)/0.8e1 + 0.48866620889e11/0.404890569012e12 * c(N-3) + 0.752806667e9/0.539854092016e12 * c(N-2) + 0.752806667e9/0.539854092016e12 * c(N-1) -0.3391e4/0.6692148e7 * c(N-4) - 0.238797444493e12/0.119037827289528e15 * c(N-3) + 0.33235054191e11/0.26452850508784e14 * c(N-2) + 0.33235054191e11/0.26452850508784e14 * c(N-1) -0.8673e4/0.2904112e7 * c(N-2) - 0.8673e4/0.2904112e7 * c(N-1) + 0.8673e4/0.1452056e7 * c(N-3) -c(N-3)/0.392e3 + c(N-2)/0.784e3 + c(N-1)/0.784e3;...
           -c(N-6)/0.6e1 - 0.1618585929605e13/0.9919818940794e13 * c(N-3) - c(N-5)/0.2e1 - 0.1118749e7/0.2230716e7 * c(N-4) - 0.13091810925e11/0.13226425254392e14 * c(N-2) - 0.13091810925e11/0.13226425254392e14 * c(N-1) c(N-6)/0.24e2 + 0.5e1/0.6e1 * c(N-5) + 0.3896014498639e13/0.4959909470397e13 * c(N-3) + 0.8386761355510099813e19/0.128413970713633903242e21 * c(N-2) + 0.280535e6/0.371786e6 * c(N-4) + 0.3360696339136261875e19/0.171218627618178537656e21 * c(N-1) -c(N-5)/0.6e1 - 0.4959271814984644613e19/0.20965546238960637264e20 * c(N-2) - 0.375177e6/0.743572e6 * c(N-4) - 0.13425842714e11/0.33740880751e11 * c(N-3) - 0.193247108773400725e18/0.6988515412986879088e19 * c(N-1) -0.365281640980e12/0.1653303156799e13 * c(N-3) + 0.564461e6/0.4461432e7 * c(N-4) + 0.1613976761032884305e19/0.7963657098519931984e19 * c(N-2) - 0.198407225513315475e18/0.7963657098519931984e19 * c(N-1) -0.1328188692663e13/0.37594290333616e14 * c(N-2) + 0.2226377963775e13/0.37594290333616e14 * c(N-1) - 0.8673e4/0.363014e6 * c(N-3) c(N-3)/0.49e2 + 0.49579087e8/0.10149031312e11 * c(N-2) - 0.256702175e9/0.10149031312e11 * c(N-1);...
           48/(49)*(-0.368395e6/0.2230716e7 * c(N-4) + c(N-5)/0.8e1 + 0.48866620889e11/0.404890569012e12 * c(N-3) + 0.752806667e9/0.539854092016e12 * c(N-2) + 0.752806667e9/0.539854092016e12 * c(N-1)) 48/(49)*(-c(N-5)/0.6e1 - 0.4959271814984644613e19/0.20965546238960637264e20 * c(N-2) - 0.375177e6/0.743572e6 * c(N-4) - 0.13425842714e11/0.33740880751e11 * c(N-3) - 0.193247108773400725e18/0.6988515412986879088e19 * c(N-1)) 48/(49)*(c(N-5)/0.24e2 + 0.1869103e7/0.2230716e7 * c(N-4) + 0.507284006600757858213e21/0.475219048083107777984e21 * c(N-2) + 0.3e1/0.1088e4 * c(N) + 0.31688435395e11/0.67481761502e11 * c(N-3) + 0.27769176016102795561e20/0.712828572124661666976e21 * c(N-1)) 48/(49)*(-0.125059e6/0.743572e6 * c(N-4) + c(N)/0.136e3 - 0.23099342648e11/0.101222642253e12 * c(N-3) - 0.4836340090442187227e19/0.5525802884687299744e19 * c(N-2) + 0.193950157930938693e18/0.5525802884687299744e19 * c(N-1)) 48/(49)*(0.260297319232891e15/0.2556411742685888e16 * c(N-2) - 0.59e2/0.1088e4 * c(N) - 0.106641839640553e15/0.1278205871342944e16 * c(N-1) + 0.26019e5/0.726028e6 * c(N-3)) 48/(49)*(-0.1244724001e10/0.21126554976e11 * c(N-2) + 0.3e1/0.68e2 * c(N) + 0.752806667e9/0.21126554976e11 * c(N-1));...
           48/(43)*(-0.3391e4/0.6692148e7 * c(N-4) - 0.238797444493e12/0.119037827289528e15 * c(N-3) + 0.33235054191e11/0.26452850508784e14 * c(N-2) + 0.33235054191e11/0.26452850508784e14 * c(N-1)) 48/(43)*(-0.365281640980e12/0.1653303156799e13 * c(N-3) + 0.564461e6/0.4461432e7 * c(N-4) + 0.1613976761032884305e19/0.7963657098519931984e19 * c(N-2) - 0.198407225513315475e18/0.7963657098519931984e19 * c(N-1)) 48/(43)*(-0.125059e6/0.743572e6 * c(N-4) + c(N)/0.136e3 - 0.23099342648e11/0.101222642253e12 * c(N-3) - 0.4836340090442187227e19/0.5525802884687299744e19 * c(N-2) + 0.193950157930938693e18/0.5525802884687299744e19 * c(N-1)) 48/(43)*(0.564461e6/0.13384296e8 * c(N-4) + 0.470299699916357e15/0.952302618316224e15 * c(N-3) + 0.550597048646198778781e21/0.1624586048098066124736e22 * c(N-1) + c(N)/0.51e2 + 0.378288882302546512209e21/0.270764341349677687456e21 * c(N-2)) 48/(43)*(-0.59e2/0.408e3 * c(N) - 0.29294615794607e14/0.29725717938208e14 * c(N-2) - 0.2234477713167e13/0.29725717938208e14 * c(N-1) - 0.8673e4/0.363014e6 * c(N-3)) 48/(43)*(-0.59e2/0.3136e4 * c(N-3) - 0.13249937023e11/0.48148892736e11 * c(N-1) + 0.2e1/0.17e2 * c(N) + 0.2083938599e10/0.8024815456e10 * c(N-2));...
           48/(59)*(-0.8673e4/0.2904112e7 * c(N-2) - 0.8673e4/0.2904112e7 * c(N-1) + 0.8673e4/0.1452056e7 * c(N-3)) 48/(59)*(-0.1328188692663e13/0.37594290333616e14 * c(N-2) + 0.2226377963775e13/0.37594290333616e14 * c(N-1) - 0.8673e4/0.363014e6 * c(N-3)) 48/(59)*(0.260297319232891e15/0.2556411742685888e16 * c(N-2) - 0.59e2/0.1088e4 * c(N) - 0.106641839640553e15/0.1278205871342944e16 * c(N-1) + 0.26019e5/0.726028e6 * c(N-3)) 48/(59)*(-0.59e2/0.408e3 * c(N) - 0.29294615794607e14/0.29725717938208e14 * c(N-2) - 0.2234477713167e13/0.29725717938208e14 * c(N-1) - 0.8673e4/0.363014e6 * c(N-3)) 48/(59)*(0.9258282831623875e16/0.7669235228057664e16 * c(N-2) + 0.3481e4/0.3264e4 * c(N) + 0.228389721191751e15/0.1278205871342944e16 * c(N-1) + 0.8673e4/0.1452056e7 * c(N-3)) 48/(59)*(-0.6025413881e10/0.21126554976e11 * c(N-2) - 0.59e2/0.68e2 * c(N) - 0.537416663e9/0.7042184992e10 * c(N-1));...
           48/(17)*(-c(N-3)/0.392e3 + c(N-2)/0.784e3 + c(N-1)/0.784e3) 48/(17)*(c(N-3)/0.49e2 + 0.49579087e8/0.10149031312e11 * c(N-2) - 0.256702175e9/0.10149031312e11 * c(N-1)) 48/(17)*(-0.1244724001e10/0.21126554976e11 * c(N-2) + 0.3e1/0.68e2 * c(N) + 0.752806667e9/0.21126554976e11 * c(N-1)) 48/(17)*(-0.59e2/0.3136e4 * c(N-3) - 0.13249937023e11/0.48148892736e11 * c(N-1) + 0.2e1/0.17e2 * c(N) + 0.2083938599e10/0.8024815456e10 * c(N-2)) 48/(17)*(-0.6025413881e10/0.21126554976e11 * c(N-2) - 0.59e2/0.68e2 * c(N) - 0.537416663e9/0.7042184992e10 * c(N-1)) 48/(17)*(0.3e1/0.3136e4 * c(N-3) + 0.27010400129e11/0.345067064608e12 * c(N-2) + 0.234566387291e12/0.690134129216e12 * c(N-1) + 0.12e2/0.17e2 * c(N))
        ];



        R = sparse(M,1);
        R(1:24) = reshape(U(1:4,:)',24,1);
        R(25:30) = U(5,:);
        R(31) = -c(5+1)/0.6e1 + c(5)/0.8e1 + c(5+2)/0.8e1;
        R(32:37) = U(6,:);
        R(38:39) = [-c(6-1)/0.6e1 - c(6+2)/0.6e1 - c(6)/0.2e1 - c(6+1)/0.2e1;...
                                      -c(6+1)/0.6e1 + c(6)/0.8e1 + c(6+2)/0.8e1];

        R(M-38:M-37) = [-c(N-6)/0.6e1 + c(N-7)/0.8e1 + c(N-5)/0.8e1;...
                                              -c(N-7)/0.6e1 - c(N-4)/0.6e1 - c(N-6)/0.2e1 - c(N-5)/0.2e1];
        R(M-36:M-31) = L(1,:);
        R(M-30) = -c(N-5)/0.6e1 + c(N-6)/0.8e1 + c(N-4)/0.8e1;
        R(M-29:M-24) = L(2,:);
        R(M-23:M) = reshape(L(3:6,:)',24,1);

        for i=7:N-6
            R(40+(i-7)*5:44+(i-7)*5) = [
                -c(i-1)/0.6e1 + c(i-2)/0.8e1 + c(i)/0.8e1,...
                -c(i-2)/0.6e1 - c(i+1)/0.6e1 - c(i-1)/0.2e1 - c(i)/0.2e1,...
                c(i-2)/0.24e2 + 0.5e1/0.6e1 * c(i-1) + 0.5e1/0.6e1 * c(i+1) + c(i+2)/0.24e2 + 0.3e1/0.4e1 * c(i),...
                -c(i-1)/0.6e1 - c(i+2)/0.6e1 - c(i)/0.2e1 - c(i+1)/0.2e1,...
                -c(i+1)/0.6e1 + c(i)/0.8e1 + c(i+2)/0.8e1
            ];
        end

        R = R/h/h;
        D2 = -R;
        D2(1:4) = -48/17/h/h*[c(1)*(-11/6);c(1)*3;c(1)*(-3/2);c(1)*1/3] + D2(1:4);
        D2(M-3:M) = -48/17/h/h*[c(N)*1/3;c(N)*(-3/2);c(N)*3;c(N)*(-11/6)] + D2(M-3:M);



%         BS = sparse(N,N);
%         BS(1,1:4) = -c(1)*1/h*[(-11/6);3;(-3/2);1/3];
%         BS(N,N-3:N) = c(N)*1/h*[(-1/3);3/2;(-3);11/6];
%         BS = sparse(BS);

    % %%Row and column indices%%
        M = 78+(N-12)*5;
        rows = [kron([1;2;3;4],ones(6,1));...
                    5*ones(7,1);...
                    6*ones(8,1);...
                    kron((7:N-6)',ones(5,1));...
                    (N-5)*ones(8,1);...
                    (N-4)*ones(7,1);...
                    kron([N-3;N-2;N-1;N],ones(6,1))];

        cols = sparse(M,1);
        cols(1:24) = kron(ones(4,1),[1;2;3;4;5;6]);
        cols(25:39) = [(1:7)';(1:8)'];
        cols(M-23:M) = kron(ones(4,1),[N-5;N-4;N-3;N-2;N-1;N]);
        cols(M-38:M-24) = [(N-7:N)';(N-6:N)'];
        for i=7:N-6
            cols(40+(i-7)*5:44+(i-7)*5) = [i-2;i-1;i;i+1;i+2];
        end
        D2 = sparse(rows,cols,D2);
    end
    D2 = @D2_fun;
end