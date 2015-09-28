function A = diagInd(d,n,m)
    A = zeros(n,length(d));
    for i = 1:length(d)
        i0 = 1;
        j0 = d(i)+1;

        I = i0 + (0:(n-1))';
        J = j0 + (0:(n-1))';

        A(:,i) = matInd2VecInd(I,J,n);

    end
end

function I = matInd2VecInd(i,j,n)
    I = i + (j-1)*n;
end