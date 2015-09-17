% Takes a step of
%   v_tt = Dv+Ev_t+S
%
%   1/k^2 * (v_next - 2v + v_prev) = Dv  + E 1/(2k)(v_next - v_prev) + S
%
function [v_next, v] = cdiff(v, v_prev, k, D, E, S)
    %   1/k^2 * (v_next - 2v + v_prev) = Dv  + E 1/(2k)(v_next - v_prev) + S
    %       ekv to
    %   A v_next = B v + C v_prev + S
    I = speye(size(D));
    A =  1/k^2 * I - 1/(2*k)*E;
    B =  2/k^2 * I + D;
    C = -1/k^2 * I - 1/(2*k)*E;

    v_next = A\(B*v + C*v_prev + S);
end