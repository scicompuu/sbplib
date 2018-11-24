% Takes a diffOp and a cell array of boundary condition definitions.
% Each bc is a struct with the fields
%  * type     -- Type of boundary condition
%  * boundary -- Boundary identifier
%  * data     -- A function_handle for a function which provides boundary data.(see below)
% Also takes S_sign which modifies the sign of the penalty function, [-1,1]
% Returns a closure matrix and a forcing function S.
%
% The boundary data function can either be a function of time or a function of time and space coordinates.
% In the case where it only depends on time it should return the data as grid function for the boundary.
% In the case where it also takes space coordinates the number of space coordinates should match the number of dimensions of the problem domain.
% For example in the 2D case: f(t,x,y).
function [closure, S] = bcSetup(diffOp, bcs, S_sign)
    default_arg('S_sign', 1);
    assertType(bcs, 'cell');
    assert(S_sign == 1 || S_sign == -1, 'S_sign must be either 1 or -1');

    [closure, penalties] = scheme.bc.closureSetup(diffOp, bcs);
    S = scheme.bc.forcingSetup(diffOp, penalties, bcs, S_sign);
end
