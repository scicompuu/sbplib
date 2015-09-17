% Converts a semidiscretized PDE from second order in time to first order in time.
%   v_tt = Dv + S
% becomes
%   w_t = Mw + C
% where
%   w = [ v ; v_t]
% and
%   M = [0 I;
%        D 0];
%   C = [0;
%        S];
function [M,C] = tt2t(D,S)
    default_arg('S',sparse(size(D)))
    time_dependent_bc = isa(S,'function_handle');

    I = eye(size(D));
    O = zeros(size(D));

    M = [O I;
         D O];

    if ~time_dependent_bc
        C = [zeros(size(S)); S];
    else
        o = zeros(size(S(0)));
        C = @(t)([o; S(t)]);
    end
end
