function [H, HI, D1, D2, e_l, e_r, d1_l, d1_r] = d2_variable_periodic_6(m,h)
    % m = number of unique grid points, i.e. h = L/m;

    if(m<7)
        error(['Operator requires at least ' num2str(7) ' grid points']);
    end

    % Norm
    Hv = ones(m,1);
    Hv = h*Hv;
    H = spdiag(Hv, 0);
    HI = spdiag(1./Hv, 0);


    % Dummy boundary operators
    e_l = sparse(m,1);
    e_r = rot90(e_l, 2);

    d1_l = sparse(m,1);
    d1_r = -rot90(d1_l, 2);


    % D1 operator
    diags   = -3:3;
    stencil = [-1/60 9/60 -45/60 0 45/60 -9/60 1/60];
    D1 = stripeMatrixPeriodic(stencil, diags, m);
    D1 = D1/h;

    % D2 operator
    scheme_width = 7;
    scheme_radius = (scheme_width-1)/2;

    r = 1:m;
    offset = scheme_width;
    r = r + offset;

    function D2 = D2_fun(c)
        c = [c(end-scheme_width+1:end); c; c(1:scheme_width) ];

        Mm3 =  c(r-2)/0.40e2 + c(r-1)/0.40e2 - 0.11e2/0.360e3 * c(r-3) - 0.11e2/0.360e3 * c(r);
        Mm2 =  c(r-3)/0.20e2 - 0.3e1/0.10e2 * c(r-1) + c(r+1)/0.20e2 + 0.7e1/0.40e2 * c(r) + 0.7e1/0.40e2 * c(r-2);
        Mm1 = -c(r-3)/0.40e2 - 0.3e1/0.10e2 * c(r-2) - 0.3e1/0.10e2 * c(r+1) - c(r+2)/0.40e2 - 0.17e2/0.40e2 * c(r) - 0.17e2/0.40e2 * c(r-1);
        M0 =  c(r-3)/0.180e3 + c(r-2)/0.8e1 + 0.19e2/0.20e2 * c(r-1) + 0.19e2/0.20e2 * c(r+1) + c(r+2)/0.8e1 + c(r+3)/0.180e3 + 0.101e3/0.180e3 * c(r);
        Mp1 = -c(r-2)/0.40e2 - 0.3e1/0.10e2 * c(r-1) - 0.3e1/0.10e2 * c(r+2) - c(r+3)/0.40e2 - 0.17e2/0.40e2 * c(r) - 0.17e2/0.40e2 * c(r+1);
        Mp2 =  c(r-1)/0.20e2 - 0.3e1/0.10e2 * c(r+1) + c(r+3)/0.20e2 + 0.7e1/0.40e2 * c(r) + 0.7e1/0.40e2 * c(r+2);
        Mp3 =  c(r+1)/0.40e2 + c(r+2)/0.40e2 - 0.11e2/0.360e3 * c(r) - 0.11e2/0.360e3 * c(r+3);

        vals = [Mm3,Mm2,Mm1,M0,Mp1,Mp2,Mp3];
        diags = -scheme_radius : scheme_radius;
        M = spdiagsPeriodic(vals,diags);

        M=M/h;
        D2=HI*(-M );
    end
    D2 = @D2_fun;


end
