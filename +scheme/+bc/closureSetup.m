% Setup closure and penalty matrices for several boundary conditions at once.
% Each bc is a struct with the fields
%  * type     -- Type of boundary condition
%  * boundary -- Boundary identifier
%  * data     -- A function_handle for a function which provides boundary data.(see below)
% Also takes S_sign which modifies the sign of the penalty function, [-1,1]
% Returns a closure matrix and a penalty matrices for each boundary condition.
%
% The boundary data function can either be a function of time or a function of time and space coordinates.
% In the case where it only depends on time it should return the data as grid function for the boundary.
% In the case where it also takes space coordinates the number of space coordinates should match the number of dimensions of the problem domain.
% For example in the 2D case: f(t,x,y).
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
