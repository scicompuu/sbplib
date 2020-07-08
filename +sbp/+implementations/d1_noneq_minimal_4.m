function [D1,H] = d1_noneq_minimal_4(N,h)

% N: Number of grid points
if(N<6)
    error('Operator requires at least 6 grid points');
end

% BP: Number of boundary points
BP = 3;

%%%% Norm matrix %%%%%%%%
P = sparse(BP,1);
%#ok<*NASGU>
P0 =  2.6864248295847e-01;
P1 =  1.0094667153500e+00;
P2 =  9.9312068011715e-01;

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
order = 4;
d = [1/12,-2/3,0,2/3,-1/12];
d = repmat(d,N,1);
Q = spdiags(d,-order/2:order/2,N,N);

% Boundaries
Q0_0 = -5.0000000000000e-01;
Q0_1 =  6.1697245625434e-01;
Q0_2 = -1.1697245625434e-01;
Q0_3 =  0.0000000000000e+00;
Q0_4 =  0.0000000000000e+00;
Q1_0 = -6.1697245625434e-01;
Q1_1 =  0.0000000000000e+00;
Q1_2 =  7.0030578958767e-01;
Q1_3 = -8.3333333333333e-02;
Q1_4 =  0.0000000000000e+00;
Q2_0 =  1.1697245625434e-01;
Q2_1 = -7.0030578958767e-01;
Q2_2 =  0.0000000000000e+00;
Q2_3 =  6.6666666666667e-01;
Q2_4 = -8.3333333333333e-02;
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