
m = 30;
h = 1;


%% 4th order non-compatible
[H, HI, D1, D2, D3, D4, e_1, e_m, M, M4,Q, Q3, S2_1, S2_m, S3_1, S3_m, S_1, S_m] = sbp.higher4(m,h);
S1 = S_1*S_1'  + S_m*S_m';
S2 = S2_1*S2_1' + S2_m*S2_m';
S3 = S3_1*S3_1' + S3_m*S3_m';

alpha_II  = util.matrixborrow(M4, h*S2  );
alpha_III = util.matrixborrow(M4, h^3*S3);
fprintf('4th order non-compatible\n')
fprintf('alpha_II:  %.10f\n',alpha_II)
fprintf('alpha_III: %.10f\n',alpha_III)
fprintf('\n')


%% 6th order non-compatible
[H, HI, D1, D2, D3, D4, e_1, e_m, M, M4,Q, Q3, S2_1, S2_m, S3_1, S3_m, S_1, S_m] = sbp.higher6(m,h);
S1 = S_1*S_1'  + S_m*S_m';
S2 = S2_1*S2_1' + S2_m*S2_m';
S3 = S3_1*S3_1' + S3_m*S3_m';

alpha_II  = util.matrixborrow(M4, h*S2  );
alpha_III = util.matrixborrow(M4, h^3*S3);
fprintf('6th order non-compatible\n')
fprintf('alpha_II:  %.10f\n',alpha_II)
fprintf('alpha_III: %.10f\n',alpha_III)
fprintf('\n')


%% 2nd order compatible
[H, HI, D1, D4, e_1, e_m, M4, Q, S2_1, S2_m, S3_1, S3_m, S_1, S_m] = sbp.higher_compatible2(m,h);
S1 = S_1*S_1'  + S_m*S_m';
S2 = S2_1*S2_1' + S2_m*S2_m';
S3 = S3_1*S3_1' + S3_m*S3_m';

alpha_II  = util.matrixborrow(M4, h*S2  );
alpha_III = util.matrixborrow(M4, h^3*S3);
fprintf('2nd order compatible\n')
fprintf('alpha_II:  %.10f\n',alpha_II)
fprintf('alpha_III: %.10f\n',alpha_III)
fprintf('\n')


%% 4th order compatible
[H, HI, D1, D4, e_1, e_m, M4, Q, S2_1, S2_m, S3_1, S3_m, S_1, S_m] = sbp.higher_compatible4(m,h);
S1 = S_1*S_1'  + S_m*S_m';
S2 = S2_1*S2_1' + S2_m*S2_m';
S3 = S3_1*S3_1' + S3_m*S3_m';

alpha_II  = util.matrixborrow(M4, h*S2  );
alpha_III = util.matrixborrow(M4, h^3*S3);
fprintf('4th order compatible\n')
fprintf('alpha_II:  %.10f\n',alpha_II)
fprintf('alpha_III: %.10f\n',alpha_III)
fprintf('\n')

%% 6th order compatible
[H, HI, D1, D4, e_1, e_m, M4, Q, S2_1, S2_m, S3_1, S3_m, S_1, S_m] = sbp.higher_compatible6(m,h);
S1 = S_1*S_1'  + S_m*S_m';
S2 = S2_1*S2_1' + S2_m*S2_m';
S3 = S3_1*S3_1' + S3_m*S3_m';

alpha_II  = util.matrixborrow(M4, h*S2  );
alpha_III = util.matrixborrow(M4, h^3*S3);
fprintf('6th order compatible\n')
fprintf('alpha_II:  %.10f\n',alpha_II)
fprintf('alpha_III: %.10f\n',alpha_III)
fprintf('\n')3





% Ordinary

for order = [2 4 6 8 10]
    op = sbp.Ordinary(m,h, order);

    S_1 = op.boundary.S_1;
    S_m = op.boundary.S_m;

    M = op.norms.M;

    S1 = S_1*S_1'  + S_m*S_m';
    alpha  = util.matrixborrow(M, h*S1);
    fprintf('%dth order Ordinary\n', order)
    fprintf('alpha:  %.10f\n', alpha)
    fprintf('\n')
end


