function [closure, penalties] = bcClosureSetup(diffOp, bcs)
    assertType(bcs, 'cell');

    % Setup storage arrays
    closure = spzeros(size(diffOp));
    penalties = cell(1, length(bcs));

    % Collect closures and penalties
    for i = 1:length(bcs)
        [localClosure, penalties{i}] = diffOp.boundary_condition(bcs{i}.boundary, bcs{i}.type);
        closure = closure + localClosure;
    end
end
