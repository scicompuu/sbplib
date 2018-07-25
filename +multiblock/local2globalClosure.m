
% Takes the local closure for ice or water and turns it into a closure for the whole system
%   local -- The local closure
%   div   -- block matrix division for the diffOp
%   I     -- Index of blockmatrix block
function closure = local2globalClosure(local, div, I)
    closure_bm = blockmatrix.zero(div);
    closure_bm{I,I} = local;

    closure = blockmatrix.toMatrix(closure_bm);
end
