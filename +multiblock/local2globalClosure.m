% Takes the block-local closures and turns it into a global closure
%   local -- The local closure
%   div   -- block matrix division for the diffOp
%   I     -- Index of blockmatrix block
function closure = local2globalClosure(local, div, I)
    closure_bm = blockmatrix.zero(div);
    closure_bm{I,I} = local;

    closure = blockmatrix.toMatrix(closure_bm);
end
