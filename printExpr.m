function printExpr(expr)
    result = evalin('caller',expr);
    fprintf('%s => %s\n', expr, toString(result));
end