function [H, HI, D1, D2, e_l, e_r, d1_l, d1_r] = d2_variable_periodic_4(m,h)
    % m = number of unique grid points, i.e. h = L/m;

    if(m<5)
        error(['Operator requires at least ' num2str(5) ' grid points']);
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

    S = d1_l*d1_l' + d1_r*d1_r';

    % D1 operator
    stencil = [1/12 -2/3 0 2/3 -1/12];
    diags = -2:2;
    Q = stripeMatrixPeriodic(stencil, diags, m);
    D1 = HI*(Q - 1/2*e_l*e_l' + 1/2*e_r*e_r');


    scheme_width = 5;
    scheme_radius = (scheme_width-1)/2;

    r = 1:m;
    offset = scheme_width;
    r = r + offset;

    function D2 = D2_fun(c)
        c = [c(end-scheme_width+1:end); c; c(1:scheme_width) ];

        % Note: these coefficients are for -M.
        Mm2 = -1/8*c(r-2) + 1/6*c(r-1) - 1/8*c(r);
        Mm1 = 1/6 *c(r-2) + 1/2*c(r-1) + 1/2*c(r) + 1/6*c(r+1);
        M0  = -1/24*c(r-2)- 5/6*c(r-1) - 3/4*c(r) - 5/6*c(r+1) - 1/24*c(r+2);
        Mp1  = 0 * c(r-2) + 1/6*c(r-1) + 1/2*c(r) + 1/2*c(r+1) + 1/6 *c(r+2);
        Mp2  = 0 * c(r-2) + 0 * c(r-1) - 1/8*c(r) + 1/6*c(r+1) - 1/8 *c(r+2);

        vals = -[Mm2,Mm1,M0,Mp1,Mp2];
        diags = -scheme_radius : scheme_radius;
        M = spdiagsPeriodic(vals,diags);

        M=M/h;
        D2=HI*(-M );

    end
    D2 = @D2_fun;
end