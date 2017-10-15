% Takes a symfun and makes a better anonymous function
function fun = matlabFunctionSizePreserving(f)
    mf = matlabFunction(f);
    args = argnames(f);

    funStr = func2str(mf);
    for i = 1:length(args)
        funStr = [funStr sprintf(' + 0*%s', toString(args(i)))];
    end

    fun = str2func(funStr);
end
