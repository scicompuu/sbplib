% Plots all the eigen values for a sparse matrix.
function plotEig(A)
    e = eig(full(A));
    e = sort(e);
    plot(e,'.')
end