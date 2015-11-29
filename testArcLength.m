m = 201; L = 1; order = 4;
close all;
tic

g1 = @(t) 2*pi*t.*sin(2*pi*t);
g2 = @(t) 2*pi*t.*cos(2*pi*t);
g = @(t) [g1(t);g2(t)];

t = linspace(0,L,m)'; dt = L/(m-1);
ops = sbp.Ordinary(m,dt,order);
D = ops.derivatives.D1;

C = grid.Curve(g,[],[],D);

% Function
figure;
C.plot([],'bo');
hold on
plot(g1(t),g2(t),'r-');
drawnow

% Derivative
figure
C.plot_derivative([],'bo');
hold on;
plot(2*pi*sin(2*pi*t) + (2*pi)^2*t.*cos(2*pi*t),...
    2*pi*cos(2*pi*t) - (2*pi)^2*t.*sin(2*pi*t),'r-')
drawnow

% Arc length
L = C.arc_length_fun(t);
figure;
plot(t,L)
drawnow

% Stretch curve
C2 = C.stretch_parameter();
z = linspace(0,1,m);
gnew = C2.g(z);
gpnew = C2.gp(z);

% Compare stretched and unstretched curves.
figure
plot(g1(t),g2(t),'b*',gnew(1,:),gnew(2,:),'ro');

% Compare stretched and unstretched derivatives. 
figure
theta = linspace(0,2*pi,100);
plot(cos(theta),sin(theta),'-b',gpnew(1,:),gpnew(2,:),'rx');

toc




