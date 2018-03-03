function [H, HI, D1, D2, e_l, e_r, d1_l, d1_r] = d2_variable_periodic_2(m,h)
    % m = number of unique grid points, i.e. h = L/m;

    if(m<3)
        error(['Operator requires at least ' num2str(3) ' grid points']);
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
    diags   = -1:1;
    stencil = [-1/2 0 1/2];
    D1 = stripeMatrixPeriodic(stencil, diags, m);
    D1 = D1/h;

    scheme_width = 3;
    scheme_radius = (scheme_width-1)/2;
    
    r = 1:m;
    offset = scheme_width;
    r = r + offset;

    function D2 = D2_fun(c)
        c = [c(end-scheme_width+1:end); c; c(1:scheme_width) ];

        Mm1 = -c(r-1)/2 - c(r)/2;
        M0  =  c(r-1)/2 + c(r)   + c(r+1)/2;
        Mp1 =            -c(r)/2 - c(r+1)/2;

        vals = [Mm1,M0,Mp1];
        diags = -scheme_radius : scheme_radius;
        M = spdiagsVariablePeriodic(vals,diags); 

        M=M/h;
        D2=HI*(-M );
    end
    D2 = @D2_fun;
end