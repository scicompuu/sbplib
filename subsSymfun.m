% Subs for a symfun
% f remains a symbolic function. If any of it's arguments is eliminated
% it is removed from the argument list while preserving the order of the
% other arguments
function f = subsSymfun(f, old, new)
    args = argnames(f);

    newExpr = subs(f, old, new);
    vars = symvar(subs(args, old, new));

    newArgs = args(ismember(args,vars));

    f = symfun(newExpr, newArgs);
end