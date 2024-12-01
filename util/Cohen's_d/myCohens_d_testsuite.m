function [new_cohensd_matrix_all] = myCohens_d_testsuite(AUC_Bench_App_TS20, ...
    AUC_Bench_App_TS40, AUC_Bench_App_TS60, AUC_Bench_App_TS80, AUC_Bench_App_TS100, bestIdx)
% myCohens_d_spsmetric function returns statistical comparison between 12 different sizes of test suites by Cohen's d effect size.
%
% Inputs:
%   AUC _Bench_App_TS20: AUC of 12 correct benchmarks under 48 different hyper configs and 7 suspiciousness metrics based on 20% test suite
%   AUC _Bench_App_TS40: AUC of 12 correct benchmarks under 48 different hyper configs and 7 suspiciousness metrics based on 40% test suite
%   AUC _Bench_App_TS60: AUC of 12 correct benchmarks under 48 different hyper configs and 7 suspiciousness metrics based on 60% test suite
%   AUC _Bench_App_TS80: AUC of 12 correct benchmarks under 48 different hyper configs and 7 suspiciousness metrics based on 80% test suite
%   AUC _Bench_App_TS100: AUC of 12 correct benchmarks under 48 different hyper configs and 7 suspiciousness metrics based on 100% test suite
%   bestIdx: 
% Outputs:
%   new_cohensd_matrix_all: statistical comparison between 12 different sizes of test suites

if ~isempty(bestIdx)
    AUC_Bench_App_TS20 = AUC_Bench_App_TS20(:,bestIdx);
    AUC_Bench_App_TS40 = AUC_Bench_App_TS40(:,bestIdx);
    AUC_Bench_App_TS60 = AUC_Bench_App_TS60(:,bestIdx);
    AUC_Bench_App_TS80 = AUC_Bench_App_TS80(:,bestIdx);
    AUC_Bench_App_TS100 = AUC_Bench_App_TS100(:,bestIdx);
end

col_AUC_Bench_App_TS20 = AUC_Bench_App_TS20(:);
col_AUC_Bench_App_TS40 = AUC_Bench_App_TS40(:);
col_AUC_Bench_App_TS60 = AUC_Bench_App_TS60(:);
col_AUC_Bench_App_TS80 = AUC_Bench_App_TS80(:);
col_AUC_Bench_App_TS100 = AUC_Bench_App_TS100(:);

% there are 5 approaches
app_num = 5;
%% consider all hyper params
AUC_Bench_App_All = [col_AUC_Bench_App_TS20, col_AUC_Bench_App_TS40, col_AUC_Bench_App_TS60, ...
    col_AUC_Bench_App_TS80, col_AUC_Bench_App_TS100];

matrix_d_All = zeros(app_num, app_num);

% calculate the Cohen's d effect size to compare each pair of apps
for i = 1:app_num
    for j = 1:app_num
        % calculate the Cohen's d effect size
        matrix_d_All(i, j) = computeCohens_d(AUC_Bench_App_All(:, i), AUC_Bench_App_All(:, j), 'paired');
    end
end

new_cohensd_matrix_all = string(matrix_d_All);
[large_worse, medium_worse, small_worse, small_better, medium_better, large_better] = label(matrix_d_All);

new_cohensd_matrix_all(large_worse) = 'large_worse';
new_cohensd_matrix_all(medium_worse) = 'medium_worse';
new_cohensd_matrix_all(small_worse) = 'small_worse';
new_cohensd_matrix_all(small_better) = 'small_better';
new_cohensd_matrix_all(medium_better) = 'medium_better';
new_cohensd_matrix_all(large_better) = 'large_better';
    function [large_worse, medium_worse, small_worse, small_better, medium_better, large_better] = label(matrix_d)
        large_worse = matrix_d < -0.8;
        medium_worse = matrix_d >= -0.8 & matrix_d <= -0.2;
        small_worse = matrix_d > -0.2 & matrix_d < 0;
        small_better = matrix_d > 0 & matrix_d < 0.2;
        medium_better = matrix_d >= 0.2 & matrix_d <= 0.8;
        large_better = matrix_d > 0.8;
    end
end