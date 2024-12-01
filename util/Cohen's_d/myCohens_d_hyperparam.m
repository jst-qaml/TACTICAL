function [new_cohensd_matrix_all] = myCohens_d_hyperparam(AUC_Bench_App_Tarantula, AUC_Bench_App_Ochiai, AUC_Bench_App_D2, AUC_Bench_App_D3, AUC_Bench_App_Jaccard, AUC_Bench_App_Kulczynski2, AUC_Bench_App_Op2, bestIdx)
% myCohens_d_hyperparam function returns statistical comparison between 48 hyper configs by Cohen's d effect size.
%
% Inputs:
%   AUC_Bench_App_Tarantula: AUC of 12 correct benchmarks under 48 different hyper configs
%   AUC_Bench_App_Ochiai:
%   AUC_Bench_App_D2:
%   AUC_Bench_App_D3:
%   AUC_Bench_App_Jaccard:
%   AUC_Bench_App_Kulczynski2:
%   AUC_Bench_App_Op2:
%   bestIdx: 
% Outputs:
%   new_cohensd_matrix_all: statistical comparison between 48 different hyper configs

if ~isempty(bestIdx)
    AUC_Bench_App_Tarantula = AUC_Bench_App_Tarantula(:,bestIdx);
    AUC_Bench_App_Ochiai = AUC_Bench_App_Ochiai(:,bestIdx);
    AUC_Bench_App_D2 = AUC_Bench_App_D2(:,bestIdx);
    AUC_Bench_App_D3 = AUC_Bench_App_D3(:,bestIdx);
    AUC_Bench_App_Jaccard = AUC_Bench_App_Jaccard(:,bestIdx);
    AUC_Bench_App_Kulczynski2 = AUC_Bench_App_Kulczynski2(:,bestIdx);
    AUC_Bench_App_Op2 = AUC_Bench_App_Op2(:,bestIdx);
end

app_num = size(AUC_Bench_App_Tarantula, 2);
%% consider all suspiciousness metrics
AUC_Bench_App_All = [AUC_Bench_App_Tarantula; AUC_Bench_App_Ochiai; AUC_Bench_App_D2; AUC_Bench_App_D3; ...
    AUC_Bench_App_Jaccard; AUC_Bench_App_Kulczynski2; AUC_Bench_App_Op2];

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