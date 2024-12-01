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
benchMutList = {
    'ACC3_spec1_simlog', ...
    'ACC3_spec2_simlog', ...
    'ACC2_spec1_simlog', ...
    'ACC2_spec2_simlog', ...
    'AFC1_spec1_simlog', ...
    'AFC1_spec2_simlog', ...
    'AFC2_spec1_simlog', ...
    'AFC2_spec2_simlog', ...
    'WT1_spec1_simlog', ...
    'WT2_spec1_simlog', ...
    'SC1_spec1_simlog', ...
    'SC2_spec1_simlog', ...
    };

nnStruList = {
    [10,10,10,10], ... % ACC3
    [10,10,10,10], ...
    [15,15,15], ...    % ACC2
    [15,15,15], ...
    [15,15,15], ...    % AFC1
    [15,15,15], ...
    [15,15,15,15], ... % AFC2
    [15,15,15,15], ...
    [5,5,5], ...       % WT1
    [15,15,15], ...    % WT2
    [10,10,10,10], ... % SC1
    [15,15,15,15], ... % SC2
    };

spsMetricList = {'Tarantula','Ochiai','D2','D3','Jaccard','Kulczynski2','Op2'};
rqIdx = 5;
% 20% faulty benchmarks for RQ1; 80% faulty benchmarks for RQ234 and RQ5.
mutPercent = 80;
% 100% test cases for RQ1234; 20%,40%,60%,80%,100% test cases for RQ5
tsPerList = [20,40,60,80,100];

for tspi = 1:5
    tsPercent = tsPerList(1,tspi);
%% count the number of faulty benchmarks
mutNumList = zeros(3,12);
for xi = 1:3
    for bmi = 1:12
        if rqIdx ~= 5
            if mutPercent == 20
                benchPath = [path, '/MUT', num2str(xi), '/', benchMutList{1,bmi}, '/RQ1'];
            elseif mutPercent == 80
                benchPath = [path, '/MUT', num2str(xi), '/', benchMutList{1,bmi}, '/RQ234'];
            end
            dirOutput = dir(fullfile(benchPath, '**', '*_M_*.mat'));
            filePaths = {dirOutput.folder};
            fileNames = {dirOutput.name};
        elseif rqIdx == 5
            benchPath = [path, '/MUT', num2str(xi), '/', benchMutList{1,bmi}, '/RQ234'];
            dirOutput = dir(fullfile(benchPath, '*_M_*.mat'));
            filePaths = {dirOutput.folder};
            fileNames = {dirOutput.name};
        end

        mutNum = numel(fileNames);
        mutNumList(xi,bmi) = mutNum;
        nnStru = nnStruList{1,bmi};
        hLNum = numel(nnStru);
        % initialize mutFLInfo and template for mutFLInfo and hyCon.
        [mutFLInfo, ~, tempHyCon] = initMutFLInfo(benchMutList{1,bmi},nnStru,mutNum,spsMetricList);
        parfor muti = 1:mutNum
            hyCon = tempHyCon;
            curMut = load(fullfile(filePaths{1,muti}, fileNames{1,muti}));
            % minOut and maxOut of curMut
            minOut = inf;
            maxOut = -inf;
            for tsi = 1:numel(curMut.cur_diagInfo_suite) * tsPercent/100
                if curMut.cur_diagInfo_suite{1,tsi}.tau_s < 1
                    continue;
                end
                matColNeuronOut = cell2mat(curMut.cur_diagInfo_suite{1,tsi}.colNeuronOut);
                tempMinOut = min(matColNeuronOut(:));
                tempMaxOut = max(matColNeuronOut(:));
                if minOut > tempMinOut
                    minOut = tempMinOut;
                end
                if maxOut < tempMaxOut
                    maxOut = tempMaxOut;
                end
            end

            for tsi = 1:numel(curMut.cur_diagInfo_suite) * tsPercent/100
                if curMut.cur_diagInfo_suite{1,tsi}.tau_s < 1
                    continue;
                end
                colNeuronOut = curMut.cur_diagInfo_suite{1,tsi}.colNeuronOut;
                state = curMut.cur_diagInfo_suite{1,tsi}.state;
                hyCon = isActivated(colNeuronOut, hyCon, nnStru, state, minOut, maxOut);
            end
            % initialize hyMCon using hyCon
            hyMCon = repmat(hyCon, numel(spsMetricList), 1);

            for mi = 1:numel(spsMetricList)
                for hci = 1:numel(hyCon)
                    % travasal the nn structure again, skip 1st hidden layer
                    idx = 0;
                    for li = 2:hLNum
                        ap = hyCon(1,hci).actSpectrum{1,li}(:,1);
                        np = hyCon(1,hci).actSpectrum{1,li}(:,2);
                        af = hyCon(1,hci).actSpectrum{1,li}(:,3);
                        nf = hyCon(1,hci).actSpectrum{1,li}(:,4);
                        spsScore = calSpsScore(ap, np, af, nf, spsMetricList{1,mi});
                        % update spsScore
                        hyMCon(mi,hci).spsNeuronList(idx+1:idx+nnStru(1,li),3) = spsScore;
                        % update idx
                        idx = idx + nnStru(1,li);
                    end
                    % I decide not to use this function
                    sortedSpsNeuronList = sortSpsNeuronList(hyMCon(mi,hci).spsNeuronList, 'fixed');
                    hyMCon(mi,hci).spsNeuronList = sortedSpsNeuronList;
                    hyMCon(mi,hci).recallList = genRecallList(sortedSpsNeuronList(:,1:2), curMut.mut_info.mut_op(:,1:2));
                end
            end
            mutFLInfo{1,muti}.hyMCon = hyMCon;
            mutFLInfo{1,muti}.mutInfo = curMut.mut_info;
            mutFLInfo{1,muti}.mutName = fileNames{1,muti};
        end

        if rqIdx ~= 5
            mutFLResName = ['Mut',num2str(xi),'_',benchMutList{1,bmi}(1:end-7),'_FL.mat'];
            if mutPercent == 20
                save(fullfile([path, '/MUT', num2str(xi), '/RQ1FL/', mutFLResName]),'mutFLInfo');
            elseif mutPercent == 80
                save(fullfile([path, '/MUT', num2str(xi), '/RQ234FL/', mutFLResName]),'mutFLInfo');
            end
        elseif rqIdx == 5
            mutFLResName = ['Mut',num2str(xi),'_',benchMutList{1,bmi}(1:end-7), '_', num2str(tsPercent), '_FL.mat'];
            save(fullfile([path, '/MUT', num2str(xi), '/RQ5FL/', mutFLResName]),'mutFLInfo');
        end
    end
end
end