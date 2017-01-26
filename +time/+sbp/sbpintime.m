% Takes one time step of size k using the sbp in time method
% starting from v_0 and where
% M v_t = E_r'*v + f(t)
function v = sbpintime(v, t, tvec, penalty, f, Nblock, E_r,...
                     L, U, P, Q)

% Pick out last time step
v_r = E_r'*v;  

% Form RHS
RHS = penalty*v_r; 

% Build vector of f-values and add to RHS
if(~isempty(f))
    fvec = [];
    for i = 1:Nblock
        fvec = [fvec;f(t+tvec(i))];
    end
    RHS = RHS + fvec;
end

% Solve system
%v = M\RHS;
RHS = P*RHS;
RHS = L\RHS;
v = U\RHS;
v = Q*v; 

