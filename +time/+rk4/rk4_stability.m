function rk_stability()
    ruku4 = @(z)(abs(1 + z +(1/2)*z.^2 + (1/6)*z.^3 + (1/24)*z.^4));
    circ  = @(z)(abs(z));


    % contour(X,Y,z)
    ax = [-4 2 -3 3];
    % hold on
    fcontour(ruku4,[1,1],[-3, 0.6],[-3.2, 3.2])
    hold on
    r = 2.6;
    fcontour(circ,[r,r],[-3, 0.6],[-3.2, 3.2],'r')
    hold off
    % contour(X,Y,z,[1,1],'b')
    axis(ax)
    title('4th order Runge-Kutta stability region')
    xlabel('Re')
    ylabel('Im')
    axis equal
    grid on
    box on
    hold off
    % surf(X,Y,z)


    rk4roots()
end

function fcontour(f,levels,x_lim,y_lim,opt)
    default_arg('opt','b')
    x = linspace(x_lim(1),x_lim(2));
    y = linspace(y_lim(1),y_lim(2));
    [X,Y] = meshgrid(x,y);
    mu = X+ 1i*Y;

    z = f(mu);

    contour(X,Y,z,levels,opt)

end


function rk4roots()
    ruku4 = @(z)(abs(1 + z +(1/2)*z.^2 + (1/6)*z.^3 + (1/24)*z.^4));
    % Roots for real evalues:
    F = @(x)(abs(ruku4(x))-1);
    real_x = fzero(F,-3);

    % Roots for imaginary evalues:
    F = @(x)(abs(ruku4(1i*x))-1);
    imag_x1 = fzero(F,-3);
    imag_x2 = fzero(F,3);


    fprintf('Real x = %f\n',real_x)
    fprintf('Imag x = %f\n',imag_x1)
    fprintf('Imag x = %f\n',imag_x2)
end