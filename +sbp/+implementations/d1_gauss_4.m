function [D1,H,x,h,e_l,e_r] = d1_gauss_4(L)

% L: Domain length
% N: Number of grid points
if(nargin < 2)
    L = 1;
end

N = 4;

% Quadrature nodes on interval [-1, 1]
x = [ -0.8611363115940526; -0.3399810435848563; 0.3399810435848563; 0.8611363115940526];

% Shift nodes to [0,L]
x = (x+1)/2*L;

% Boundary extrapolation operators
e_l = [1.5267881254572668; -0.8136324494869273; 0.4007615203116504; -0.1139171962819899];
e_r = flipud(e_l);
e_l = sparse(e_l);
e_r = sparse(e_r);

%%%% Compute approximate h %%%%%%%%%%
h = L/(N-1);
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Norm matrix on [-1,1] %%%%%%%%
P = sparse(N,N);
P(1,1) =  0.3478548451374539;
P(2,2) =  0.6521451548625461;
P(3,3) =  0.6521451548625461;
P(4,4) =  0.3478548451374539;
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Norm matrix on [0,L] %%%%%%%%
H = P*L/2;
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% D1 on [-1,1] %%%%%%%%
D1 = sparse(N,N);
D1(1,1) = -3.3320002363522817;
D1(1,2) = 4.8601544156851962;
D1(1,3) = -2.1087823484951789;
D1(1,4) = 0.5806281691622644;

D1(2,1) = -0.7575576147992339;
D1(2,2) = -0.3844143922232086;
D1(2,3) = 1.4706702312807167;
D1(2,4) = -0.3286982242582743;

D1(3,1) = 0.3286982242582743;
D1(3,2) = -1.4706702312807167;
D1(3,3) = 0.3844143922232086;
D1(3,4) = 0.7575576147992339;

D1(4,1) = -0.5806281691622644;
D1(4,2) = 2.1087823484951789;
D1(4,3) = -4.8601544156851962;
D1(4,4) = 3.3320002363522817;
%%%%%%%%%%%%%%%%%%%%%%%%%

% D1 on [0,L]
D1 = D1*2/L;