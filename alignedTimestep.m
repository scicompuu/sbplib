% Calcualtes the largest timestep smaller than k_max that gives an integer
% number of timesteps to time T.
%   k_max -- largest allowable timestep
%   T     -- time to align with
%
% Returns:
%   k -- calculated timestep
%   N -- number to of timestep to reach T
function [k, N] = alignedTimestep(k_max, T)
    N = ceil(T/k_max);
    k = T/N;
end