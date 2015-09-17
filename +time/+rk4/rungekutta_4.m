% Takes one time step of size k using the rungekutta method
% starting from v_0 and where the function F(v,t) gives the
% time derivatives.
function v = rungekutta_4(v, t , k, F)
    k1 = F(v         ,t      );
    k2 = F(v+0.5*k*k1,t+0.5*k);
    k3 = F(v+0.5*k*k2,t+0.5*k);
    k4 = F(v+    k*k3,t+    k);
    v = v + (1/6)*(k1+2*(k2+k3)+k4)*k;
end