function axPos(n,m)
    a = {};
    aPos = [];
    for i = 1:n*m
        a{i} = subplot(n, m, i);
        aPos(i,:) = a{i}.Position;
    end

    dx = aPos(1,3);
    dy = aPos(1,4);

    x = unique(aPos(:,1));
    y = unique(aPos(:,2));


    fprintf('dx: %f\n', dx);
    fprintf('dy: %f\n', dy);

    fprintf(' x: ');
    fprintf('%f ', x);
    fprintf('\n');

    fprintf(' y: ');
    fprintf('%f ', y);
    fprintf('\n');
end
