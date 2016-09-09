function [D1,H,x,h] = d1_noneq_4(N,L)

% L: Domain length
% N: Number of grid points
if(nargin < 2)
    L = 1;
end

if(N<8)
    error('Operator requires at least 8 grid points');
end

% BP: Number of boundary points
% m:  Number of nonequidistant spacings
% order: Accuracy of interior stencil
BP = 4;
m = 2;
order = 4;

%%%% Non-equidistant grid points %%%%%
x0 =  0.0000000000000e+00;
x1 =  6.8764546205559e-01;
x2 =  1.8022115125776e+00;
x3 =  2.8022115125776e+00;
x4 =  3.8022115125776e+00;

xb = zeros(m+1,1);
for i = 0:m
    xb(i+1) = eval(['x' num2str(i)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Compute h %%%%%%%%%%
h = L/(2*xb(end) + N-1-2*m);
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Define grid %%%%%%%%
x = h*[xb; linspace(xb(end)+1,L/h-xb(end)-1,N-2*(m+1))'; L/h-flip(xb) ];
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Norm matrix %%%%%%%%
P = zeros(BP,1);
%#ok<*NASGU>
P0 =  2.1259737557798e-01;
P1 =  1.0260290400758e+00;
P2 =  1.0775123588954e+00;
P3 =  9.8607273802835e-01;

for i = 0:BP-1
    P(i+1) = eval(['P' num2str(i)]);
end

H = ones(N,1);
H(1:BP) = P;
H(end-BP+1:end) = flip(P);
H = spdiags(h*H,0,N,N);
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Q matrix %%%%%%%%%%%

% interior stencil
switch order
    case 2
        d = [-1/2,0,1/2];
    case 4
        d = [1/12,-2/3,0,2/3,-1/12];
    case 6
        d = [-1/60,3/20,-3/4,0,3/4,-3/20,1/60];
    case 8
        d = [1/280,-4/105,1/5,-4/5,0,4/5,-1/5,4/105,-1/280];
    case 10
        d = [-1/1260,5/504,-5/84,5/21,-5/6,0,5/6,-5/21,5/84,-5/504,1/1260];
    case 12
        d = [1/5544,-1/385,1/56,-5/63,15/56,-6/7,0,6/7,-15/56,5/63,-1/56,1/385,-1/5544];
end
d = repmat(d,N,1);
Q = spdiags(d,-order/2:order/2,N,N);

% Boundaries
Q0_0 = -5.0000000000000e-01;
Q0_1 =  6.5605279837843e-01;
Q0_2 = -1.9875859409017e-01;
Q0_3 =  4.2705795711740e-02;
Q0_4 =  0.0000000000000e+00;
Q0_5 =  0.0000000000000e+00;
Q1_0 = -6.5605279837843e-01;
Q1_1 =  0.0000000000000e+00;
Q1_2 =  8.1236966439895e-01;
Q1_3 = -1.5631686602052e-01;
Q1_4 =  0.0000000000000e+00;
Q1_5 =  0.0000000000000e+00;
Q2_0 =  1.9875859409017e-01;
Q2_1 = -8.1236966439895e-01;
Q2_2 =  0.0000000000000e+00;
Q2_3 =  6.9694440364211e-01;
Q2_4 = -8.3333333333333e-02;
Q2_5 =  0.0000000000000e+00;
Q3_0 = -4.2705795711740e-02;
Q3_1 =  1.5631686602052e-01;
Q3_2 = -6.9694440364211e-01;
Q3_3 =  0.0000000000000e+00;
Q3_4 =  6.6666666666667e-01;
Q3_5 = -8.3333333333333e-02;
for i = 1:BP
    for j = 1:BP
        Q(i,j) = eval(['Q' num2str(i-1) '_' num2str(j-1)]);
        Q(N+1-i,N+1-j) = -eval(['Q' num2str(i-1) '_' num2str(j-1)]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Difference operator %%
D1 = H\Q;
%%%%%%%%%%%%%%%%%%%%%%%%%%%