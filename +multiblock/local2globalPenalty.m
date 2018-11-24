% Takes the block-local penalty and turns it into a global penalty
%   local -- The local penalty
%   div   -- block matrix division for the diffOp
%   I     -- Index of blockmatrix block
function penalty = local2globalPenalty(local, div, I)
    penaltyDiv = {div{1}, size(local,2)};
    penalty_bm = blockmatrix.zero(penaltyDiv);
    penalty_bm{I,1} = local;

    penalty = blockmatrix.toMatrix(penalty_bm);
end
