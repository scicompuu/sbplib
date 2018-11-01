% Setup closure and penalty matrices for several boundary conditions at once.
function [closure, penalties] = closureSetup(diffOp, bcs)
    scheme.bc.verifyFormat(bcs, diffOp);

    % Setup storage arrays
    closure = spzeros(size(diffOp));
    penalties = cell(1, length(bcs));

    % Collect closures and penalties
    for i = 1:length(bcs)
        [localClosure, penalties{i}] = diffOp.boundary_condition(bcs{i}.boundary, bcs{i}.type);
        closure = closure + localClosure;
    end
end
