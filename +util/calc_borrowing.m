function calc_borrowing(m, h)
    default_arg('m',100);
    default_arg('h',1);

    operators = {
        {
            'd4_lonely', getM4_lonely, {
                {4, 'min_boundary_points'},
                {6, 'min_boundary_points'},
                {6, '2'},
                {6, '3'},
                {8, 'min_boundary_points'},
                {8, 'higher_boundary_order'},
            }
        }, {
            'd4_variable', {
                {2},
                {4},
                {6},
            }
        }
        % BORKEN BAD IDEA
    }


    for i = 1:operators
        baseName = operators{i}{1};
        postFixes = operators{i}{2};
        for pf = postFixes
            [a2, a3] = borrowFromD4(m, h, l{:});
        end
    end



    lonely = {
        {4, 'min_boundary_points'},
        {6, 'min_boundary_points'},
        {6, '2'},
        {6, '3'},
        {8, 'min_boundary_points'},
        {8, 'higher_boundary_order'},
    };

    for i = 1:length(lonely)
        l = lonely{i};
        [a2, a3] = d4_lonely(m, h, l{:});
        fprintf('d4_lonely %d %s\n', l{:})
        fprintf('\t  alpha_II = %f\n', a2)
        fprintf('\t alpha_III = %f\n', a3)
        fprintf('\n')
    end

    variable = {
        {2},
        {4},
        {6},
    };

    for i = 1:length(variable)
        l = variable{i};
        [a2, a3] = d4_variable(m, h, l{:});
        fprintf('d4_variable %d\n', l{:})
        fprintf('\t  alpha_II = %f\n', a2)
        fprintf('\t alpha_III = %f\n', a3)
        fprintf('\n')
    end


    %% 4th order non-compatible
    [H, HI, D1, D2, D3, D4, e_1, e_m, M, M4,Q, Q3, S2_1, S2_m, S3_1, S3_m, S_1, S_m] = sbp.higher4(m,h);
    S1 = S_1*S_1'  + S_m*S_m';
    S2 = S2_1*S2_1' + S2_m*S2_m';
    S3 = S3_1*S3_1' + S3_m*S3_m';

    alpha_I  = util.matrixborrow(M4, h^-1*S1  );
    alpha_II  = util.matrixborrow(M4, h*S2  );
    alpha_III = util.matrixborrow(M4, h^3*S3);
    fprintf('4th order non-compatible\n')
    fprintf('alpha_I1:  %.10f\n',alpha_I)
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
    fprintf('\n')





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




end

function [alpha_II, alpha_III] = d4_lonely(m, h, order, modifier)
    default_arg('modifier', [])
    func = sprintf('sbp.implementations.d4_lonely_%d', order);
    if ~isempty(modifier)
        func = sprintf('%s_%s', func, modifier);
    end
    funcCall = sprintf('%s(%s,%s)', func, toString(m), toString(h));
    [H, HI, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = eval(funcCall);

    d2d2 = d2_l*d2_l' + d2_r*d2_r';
    alpha_II  = util.matrixborrow(M4, h*d2d2);

    d3d3 = d3_l*d3_l' + d3_r*d3_r';
    alpha_III = util.matrixborrow(M4, h^3*d3d3);
end

function [alpha_II, alpha_III] = d4_variable(m, h, order)
    default_arg('modifier', [])
    func = sprintf('sbp.implementations.d4_variable_%d', order);

    funcCall = sprintf('%s(%s,%s)', func, toString(m), toString(h));
    [H, HI, D1, D2, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = eval(funcCall);

    d2d2 = d2_l*d2_l' + d2_r*d2_r';
    alpha_II  = util.matrixborrow(M4, h*d2d2);

    d3d3 = d3_l*d3_l' + d3_r*d3_r';
    alpha_III = util.matrixborrow(M4, h^3*d3d3);
end

function [d2_l, d2_r, d3_l, d3_r, M4] = getM4_lonely(m, h, order, modifier)
    fStr = getFunctionCallStr('d4_lonely', {order, modifier}, {m ,h});
    [H, HI, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = eval(funcCall);
end


% Calculates the borrowing constants for a D4 operator.
% getM4 is a function handle on the form
%  [d2_l, d2_r, d3_l, d3_r, M4] = getM4(m,h)
function [a2, a3] = borrowFromD4(m, h, getM4)
    [d2_l, d2_r, d3_l, d3_r, M4] = getM4(m, h);

    d2d2 = d2_l*d2_l' + d2_r*d2_r';
    a2  = util.matrixborrow(M4, h*d2d2);

    d3d3 = d3_l*d3_l' + d3_r*d3_r';
    a3 = util.matrixborrow(M4, h^3*d3d3);
end


function funcCallStr = getFunctionCallStr(baseName, postFix, parameters)
    default_arg('postFix', [])
    default_arg('parameters', [])

    funcCallStr = sprintf('sbp.implementations.%s', baseName);

    for i = 1:length(postFix)
        if ischar(postFix{i})
            funcCallStr = [funcCallStr '_' postFix{i}];
        else
            funcCallStr = [funcCallStr '_' toString(postFix{i})];
        end
    end

    if isempty(parameters)
        return
    end

    funcCallStr = [funcCallStr '(' toString(parameters{1})];

    for i = 2:length(parameters)
        funcCallStr = [funcCallStr ', ' toString(parameters{i})];
    end

    funcCallStr = [funcCallStr ')';
end
