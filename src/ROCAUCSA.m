% This function can extract ROC and AUC from mutFLInfo in MutX_bench_FL.mat
% and perform statistical analysis for RQ12345.
% There are three hyper parameters:
% threshold(h):          {0,5%,10%}
% top(k):                for WT#1:{1,2,3}; for ACC#1/SC#1:{2,3,4}; for the rest:{3,4,5}
% time interval(deltaT): {0.5,1.0,1.5}

% INA (h) sml  {0,5%,10%}
% ITK (k) sml  {1,2,3}/{2,3,4}/{3,4,5}
% PNA (h, deltaT) ss sm sl ms mm ml ls lm ll  {0,5%,10%} * {0.5,1.0,1.5}
% PTK (deltaT, k) ss sm sl ms mm ml ls lm ll  {0.5,1.0,1.5} * {1,2,3}/{2,3,4}/{3,4,5}
% PD  (h, deltaT) ss sm sl ms mm ml ls lm ll  {0,5%,10%} * {0.5,1.0,1.5}
% ND  (h, deltaT) ss sm sl ms mm ml ls lm ll  {0,5%,10%} * {0.5,1.0,1.5}
% MI  (deltaT) sml  {0.5,1.0,1.5}
% MD  (deltaT) sml  {0.5,1.0,1.5}
clear;
clc;
close all;
%%
path = '/Users/ldy/git/TACTICAL';
addpath(genpath(path));
benchMutFLInfoList = {
    'MutX_ACC3_spec1_FL', ...
    'MutX_ACC3_spec2_FL', ...
    'MutX_ACC2_spec1_FL', ...
    'MutX_ACC2_spec2_FL', ...
    'MutX_AFC1_spec1_FL', ...
    'MutX_AFC1_spec2_FL', ...
    'MutX_AFC2_spec1_FL', ...
    'MutX_AFC2_spec2_FL', ...
    'MutX_WT1_spec1_FL', ...
    'MutX_WT2_spec1_FL', ...
    'MutX_SC1_spec1_FL', ...
    'MutX_SC2_spec1_FL', ...
    };

topsList = [30,30,...
    30,30, ...
    30,30, ...
    45,45, ...
    10, ...
    30, ...
    30, ...
    45];

spsMetricList = {'Tarantula','Ochiai','D2','D3','Jaccard','Kulczynski2','Op2'};
metricNum = numel(spsMetricList);

% rqIdx-mutPercent-tsPercent param pairs:
% RQ1: <1,20>; RQ23: <23,80>; RQ4: <4,80>; RQ5: <5,80,20/40/60/80/100>
rqIdx = 1;
if rqIdx == 1
    mutPercent = 20;
elseif rqIdx == 23 || rqIdx == 4 
    mutPercent = 80;
elseif rqIdx == 5
    mutPercent = 80;
    % 20/40/60/80/100
    tsPercent = 20;
end

benchHyMConRecallList = cell(12,1);
mutNumList = zeros(3,12);

for bmi = 1:12
    % store recallList of MUTx (x = 1,2,3)
    hyMConRecallList = cell(metricNum, 48);
    % aggregate MUTx (x = 1,2,3)
    for xi = 1:3
        if rqIdx ~= 5
            if mutPercent == 20
                mutFLInfoFile = [path, '/MUT', num2str(xi), '/RQ1FL/', strrep(benchMutFLInfoList{1,bmi}, 'X', num2str(xi)), '.mat'];
            elseif mutPercent == 80
                mutFLInfoFile = [path, '/MUT', num2str(xi), '/RQ234FL/', strrep(benchMutFLInfoList{1,bmi}, 'X', num2str(xi)), '.mat'];
            end
        elseif rqIdx == 5 && mutPercent == 80
            mutFLInfoPrefix = strrep(benchMutFLInfoList{1,bmi}, 'X', num2str(xi));
            mutFLInfoFile = [path, '/MUT', num2str(xi), '/RQ5FL/', insertAfter(mutFLInfoPrefix, strlength(mutFLInfoPrefix) - 3, ['_', num2str(tsPercent)]), '.mat'];
        end

        % load mutFLInfo
        load(mutFLInfoFile);
        mutNum = numel(mutFLInfo);
        mutNumList(xi, bmi) = mutNum;
        % for current xi
        for muti = 1:mutNum
            for mi = 1:metricNum
                for hci = 1:48
                    % recallList of curMut as a new column
                    hyMConRecallList{mi,hci}(:,end+1) = mutFLInfo{1,muti}.hyMCon(mi,hci).recallList;
                end
            end
        end
    end
    benchHyMConRecallList{bmi,1} = hyMConRecallList;
end

% aggregate detection rate of each mutant
benchHyMConRecall = cell(12,1);
for bmi = 1:12
    hyMConRecallList = benchHyMConRecallList{bmi,1};
    hyMConRecall = cell(metricNum,48);
    maxTops = topsList(1,bmi);
    for mi = 1:metricNum
        for hci = 1:48
            detectRate = zeros(maxTops,1);
            for topsi = 1:maxTops
                % if topsi == 1
                %     detectRate(topsi,1) = mean(hyMConRecallList{mi,hci}(1,1:mutNumList(1,bmi)));
                % elseif topsi == 2
                %     detectRate(topsi,1) = mean(hyMConRecallList{mi,hci}(2,1:mutNumList(1,bmi)+mutNumList(2,bmi)));
                % else
                detectRate(topsi,1) = mean(hyMConRecallList{mi,hci}(topsi,:));
                % end
            end
            hyMConRecall{mi,hci} = detectRate;
        end
    end
    benchHyMConRecall{bmi,1} = hyMConRecall;
end

% calculate AUC, normalize Tops
benchHyMConAUC = cell(12,1);
for bmi = 1:12
    hyMConAUC = zeros(metricNum,48);
    for mi = 1:metricNum
        for hci = 1:48
            maxTops = topsList(1,bmi);
            curRecall = benchHyMConRecall{bmi,1}{mi,hci};
            hyMConAUC(mi,hci) = 5 * trapz([0, 1:1:maxTops/5], [0; curRecall(1:maxTops/5,1)]')/maxTops;
        end
    end
    benchHyMConAUC{bmi,1} = hyMConAUC;
end

aucTarantula = zeros(12,48);
aucOchiai = zeros(12,48);
aucD2 = zeros(12,48);
aucD3 = zeros(12,48);
aucJaccard = zeros(12,48);
aucKulczynski2 = zeros(12,48);
aucOp2 = zeros(12,48);

for bmi = 1:12
    aucTarantula(bmi,:) = benchHyMConAUC{bmi,1}(1,:);
    aucOchiai(bmi,:) = benchHyMConAUC{bmi,1}(2,:);
    aucD2(bmi,:) = benchHyMConAUC{bmi,1}(3,:);
    aucD3(bmi,:) = benchHyMConAUC{bmi,1}(4,:);
    aucJaccard(bmi,:) = benchHyMConAUC{bmi,1}(5,:);
    aucKulczynski2(bmi,:) = benchHyMConAUC{bmi,1}(6,:);
    aucOp2(bmi,:) = benchHyMConAUC{bmi,1}(7,:);
end

if rqIdx == 23
    aucTarantula(:,end+1) = 0.1 * ones(12,1);
    aucOchiai(:,end+1) = 0.1 * ones(12,1);
    aucD2(:,end+1) = 0.1 * ones(12,1);
    aucD3(:,end+1) = 0.1 * ones(12,1);
    aucJaccard(:,end+1) = 0.1 * ones(12,1);
    aucKulczynski2(:,end+1) = 0.1 * ones(12,1);
    aucOp2(:,end+1) = 0.1 * ones(12,1);
end

aucAll = [aucTarantula; aucOchiai; aucD2; aucD3; aucJaccard; aucKulczynski2; aucOp2];

if rqIdx == 1
    bestCrIdx = [];
    % find the best hyper params for each criterion, 20% faulty benchmarks
    [new_cohensd_matrix_all_hyperparam] = myCohens_d_hyperparam(aucTarantula, aucOchiai, aucD2, aucD3, aucJaccard, aucKulczynski2, aucOp2, bestCrIdx);
    bestCrIdx = findBestCrIdx(new_cohensd_matrix_all_hyperparam);
elseif rqIdx == 23
    % [INA, ITK, PNA, PTK, PD, ND, MI, MD]
    bestCrIdx = [1,4,7,20,27,34,43,46,49]; % bestIdx from RQ1 (20% faulty benchmarks)
    % use the best hyper params for each criterion, 80% faulty benchmarks
    [new_cohensd_matrix_all_hyperparam] = myCohens_d_hyperparam(aucTarantula, aucOchiai, aucD2, aucD3, aucJaccard, aucKulczynski2, aucOp2, bestCrIdx);
elseif rqIdx == 4
    % [INA, ITK, PNA, PTK, PD, ND, MI, MD]
    bestCrIdx = [1,4,7,20,27,34,43,46];    % bestIdx from RQ1 (20% faulty benchmarks)
    % use the best hyper params for each criterion, 80% faulty benchmarks
    [new_cohensd_matrix_all_spsmetric] = myCohens_d_spsmetric(aucTarantula, aucOchiai, aucD2, aucD3, aucJaccard, aucKulczynski2, aucOp2, bestCrIdx);
    bestMetIdx = findBestMetIdx(new_cohensd_matrix_all_spsmetric);
elseif rqIdx == 5
    bestCrIdx = [];
    % find best hyper param for each criterion, randomly selected faulty benchmarks
    [new_cohensd_matrix_all_hyperparam] = myCohens_d_hyperparam(aucTarantula, aucOchiai, aucD2, aucD3, aucJaccard, aucKulczynski2, aucOp2, bestCrIdx);
    bestCrIdx = findBestCrIdx(new_cohensd_matrix_all_hyperparam);
    [new_cohensd_matrix_all_spsmetric] = myCohens_d_spsmetric(aucTarantula, aucOchiai, aucD2, aucD3, aucJaccard, aucKulczynski2, aucOp2, bestCrIdx);
    bestMetIdx = findBestMetIdx(new_cohensd_matrix_all_spsmetric);
end

% statistical analysis results
if rqIdx == 1 || rqIdx == 23
    saResName = ['saRes_RQ', num2str(rqIdx), '_', num2str(mutPercent), '%_bestCrIdx_', num2str(~isempty(bestCrIdx)), '.mat'];
elseif rqIdx == 4
    saResName = ['saRes_RQ', num2str(rqIdx), '_', num2str(mutPercent), '%_bestCrIdx_', num2str(~isempty(bestCrIdx)), '%_bestMetIdx_', num2str(isempty(bestMetIdx)), '.mat'];
elseif rqIdx == 5
    saResName = ['saRes_RQ', num2str(rqIdx), '_', num2str(mutPercent), '%_', num2str(tsPercent), '%_bestCrIdx_', num2str(~isempty(bestCrIdx)), '_bestMetIdx_', num2str(~isempty(bestMetIdx)), '.mat'];
end

if rqIdx == 1 || rqIdx == 23
    save(fullfile([path, '/statistical analysis/RQ', num2str(rqIdx), '/', saResName]), 'rqIdx', 'mutPercent', 'bestCrIdx', 'new_cohensd_matrix_all_hyperparam');
elseif rqIdx == 4
    save(fullfile([path, '/statistical analysis/RQ4/', saResName]), 'rqIdx', 'mutPercent', 'bestCrIdx', 'bestMetIdx', 'new_cohensd_matrix_all_spsmetric');
elseif rqIdx == 5
    save(fullfile([path, '/statistical analysis/RQ', num2str(rqIdx), '/', saResName]), 'rqIdx', 'mutPercent', 'tsPercent', 'bestCrIdx', 'bestMetIdx', 'aucAll', ...
        'benchHyMConAUC','new_cohensd_matrix_all_hyperparam');
end
%% find bestCrIdx 
function [bestCrIdx] = findBestCrIdx(new_cohensd_matrix_all_hyperparam)
bestCrIdx = zeros(1,8);
inaSE = [1,3];
itkSE = [4,6];
pnaSE = [7,15];
ptkSE = [16,24];
pdSE = [25,33];
ndSE = [34,42];
miSE = [43,45];
mdSE = [46,48];
separator = [inaSE; itkSE; pnaSE; ptkSE; pdSE; ndSE; miSE; mdSE];

for ci = 1:8
    roi_new_cohensd_matrix_all_hyperparam = new_cohensd_matrix_all_hyperparam(separator(ci,1):separator(ci,2), separator(ci,1):separator(ci,2));
    for i = 1:size(roi_new_cohensd_matrix_all_hyperparam,1)
        isBest = contains(roi_new_cohensd_matrix_all_hyperparam(i,:), 'better');
        if (sum(isBest == 0) == 1)
            bestCrIdx(1,ci) = separator(ci,1) + i - 1;
            break;
        end
    end
end
end
%% find bestMetIdx 
function [bestMetIdx] = findBestMetIdx(new_cohensd_matrix_all_spsmetric)
for i = 1:size(new_cohensd_matrix_all_spsmetric,1)
    isBest = contains(new_cohensd_matrix_all_spsmetric(i,:), 'better');
    if (sum(isBest == 0) == 1)
        bestMetIdx = i;
        break;
    end
end
end
%% This section performs comparison between test suites of different sizes using statistical analysis
RQ5ResFolder = '/Users/ldy/git/TACTICAL/statistical analysis/RQ5';
tsPerList = [20,40,60,80,100];
tsPerNum = length(tsPerList);
aucAllList = cell(1,tsPerNum);
for tspi = 1:numel(tsPerList)
    resFile = ['saRes_RQ5_80%_', num2str(tsPerList(1,tspi)), '%_bestCrIdx_1_bestMetIdx_1.mat'];
    curSARes = load(resFile);
    % here, we only use Talantula
    aucAllList{1,tspi} = curSARes.aucAll(1:12,curSARes.bestCrIdx);
end

[new_cohensd_matrix_all_testsuite] = myCohens_d_testsuite(aucAllList{1,1}, aucAllList{1,2}, aucAllList{1,3}, aucAllList{1,4}, ...
    aucAllList{1,5}, []);
