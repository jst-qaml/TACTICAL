function [mutFLInfo, hyMCon, hyCon] = initMutFLInfo(bench, nnStru, mutNum, spsMetricList)
% initMutFLInfo function returns initialized mutFLInfo, hyMCon and hyCon.
%
% Inputs:
%   bench: benchmark
%   nnStru: nn structure
%   mutNum: the number of mutants
%   spsMetricList: suspiciousness metrics used for FL
% Outputs:
%   mutFLInfo: an empty fault localisation results of all mutants
%   hyMCon: an empty struct array
%   hyCon: an empty struct array

% initialize the mutFLInfo
metricNum = numel(spsMetricList);
hLNum = numel(nnStru);
% initialize the actSpectrum based on current nn structure
actSpectrum = cell(1,hLNum);
layerActSpectrum = zeros(nnStru(1,1),4);
for li = 2:hLNum
    actSpectrum{1,li} = layerActSpectrum;
end

% li,i,spsScore
spsNeuronList = [];
for li = 2:hLNum
    for j = 1:nnStru(1,li)
        spsNeuronList(end+1,:) = [li,j,0];
    end
end

% hyper config (hycon) includes 48 hyper-param combinations and their corresponding FL results
hyCon(1,48) = struct( ...
    'criterion', [], ...
    'hypername1', [], ...
    'hypername2', [], ...
    'hypersize1', [], ...
    'hypersize2', [], ...
    'hyperval1', [], ...
    'hyperval2', [], ...
    'spsNeuronList', [], ...
    'actSpectrum', [], ...
    'recallList', []);

% assign actSpectrum to each struct in hyCon
for hci = 1:48
    hyCon(hci).spsNeuronList = spsNeuronList;
    hyCon(hci).actSpectrum = actSpectrum;  % directly assign actSpectrum to each struct
end

% INA
[hyCon(1,1:3).criterion] = deal('INA');
[hyCon(1,1:3).hypername1] = deal('h');
[hyCon(1,1:3).hypersize1] = deal('S','M','L');
[hyCon(1,1:3).hyperval1] = deal(0,0.05,0.1);
% ITK
[hyCon(1,4:6).criterion] = deal('ITK');
[hyCon(1,4:6).hypername1] = deal('k');
[hyCon(1,4:6).hypersize1] = deal('S','M','L');
if contains(bench, 'WT1')
    [hyCon(1,4:6).hyperval1] = deal(1,2,3);
elseif contains(bench, 'ACC3') || contains(bench, 'SC1')
    [hyCon(1,4:6).hyperval1] = deal(2,3,4);
else
    [hyCon(1,4:6).hyperval1] = deal(3,4,5);
end
% PNA
[hyCon(1,7:15).criterion] = deal('PNA');
[hyCon(1,7:15).hypername1] = deal('h');
[hyCon(1,7:15).hypername2] = deal('deltaT');
[hyCon(1,7:15).hypersize1] = deal('S','S','S','M','M','M','L','L','L');
[hyCon(1,7:15).hypersize2] = deal('S','M','L','S','M','L','S','M','L');
[hyCon(1,7:15).hyperval1] = deal(0,0,0,  0.05,0.05,0.05,  0.1,0.1,0.1);
[hyCon(1,7:15).hyperval2] = deal(0.5,1.0,1.5,  0.5,1.0,1.5,  0.5,1.0,1.5);
% PTK
[hyCon(1,16:24).criterion] = deal('PTK');
[hyCon(1,16:24).hypername1] = deal('deltaT');
[hyCon(1,16:24).hypername2] = deal('k');
[hyCon(1,16:24).hypersize1] = deal('S','S','S','M','M','M','L','L','L');
[hyCon(1,16:24).hypersize2] = deal('S','M','L','S','M','L','S','M','L');
[hyCon(1,16:24).hyperval1] = deal(0.5,0.5,0.5,  1.0,1.0,1.0,  1.5,1.5,1.5);
if contains(bench, 'WT1')
    [hyCon(1,16:24).hyperval2] = deal(1,2,3,  1,2,3,  1,2,3);
elseif contains(bench, 'ACC3') || contains(bench, 'SC1')
    [hyCon(1,16:24).hyperval2] = deal(2,3,4,  2,3,4,  2,3,4);
else
    [hyCon(1,16:24).hyperval2] = deal(3,4,5,  3,4,5,  3,4,5);
end
% PD
[hyCon(1,25:33).criterion] = deal('PD');
[hyCon(1,25:33).hypername1] = deal('h');
[hyCon(1,25:33).hypername2] = deal('deltaT');
[hyCon(1,25:33).hypersize1] = deal('S','S','S','M','M','M','L','L','L');
[hyCon(1,25:33).hypersize2] = deal('S','M','L','S','M','L','S','M','L');
[hyCon(1,25:33).hyperval1] = deal(0,0,0,  0.05,0.05,0.05,  0.1,0.1,0.1);
[hyCon(1,25:33).hyperval2] = deal(0.5,1.0,1.5,  0.5,1.0,1.5,  0.5,1.0,1.5);
% ND
[hyCon(1,34:42).criterion] = deal('ND');
[hyCon(1,34:42).hypername1] = deal('h');
[hyCon(1,34:42).hypername2] = deal('deltaT');
[hyCon(1,34:42).hypersize1] = deal('S','S','S','M','M','M','L','L','L');
[hyCon(1,34:42).hypersize2] = deal('S','M','L','S','M','L','S','M','L');
[hyCon(1,34:42).hyperval1] = deal(0,0,0,  0.05,0.05,0.05,  0.1,0.1,0.1);
[hyCon(1,34:42).hyperval2] = deal(0.5,1.0,1.5,  0.5,1.0,1.5,  0.5,1.0,1.5);
% MI
[hyCon(1,43:45).criterion] = deal('MI');
[hyCon(1,43:45).hypername1] = deal('deltaT');
[hyCon(1,43:45).hypersize1] = deal('S','M','L');
[hyCon(1,43:45).hyperval1] = deal(0.5,1.0,1.5);
% MD
[hyCon(1,46:48).criterion] = deal('MD');
[hyCon(1,46:48).hypername1] = deal('deltaT');
[hyCon(1,46:48).hypersize1] = deal('S','M','L');
[hyCon(1,46:48).hyperval1] = deal(0.5,1.0,1.5);

hyMCon = repmat(hyCon, metricNum, 1);

% design new data structure to store the FL results of each mutant
% actually, we only need topmax, i.e., the last cell of spsneuronInfoCell.
mutFLInfo = cell(1,mutNum);
curMutInfo = struct('hyMCon', hyMCon, 'mutInfo', [], 'mutName', []);
% for muti = 1:numel(spsneuronInfoCell{1,1})/mutGen
for muti = 1:mutNum
    mutFLInfo{1,muti} = curMutInfo;
end
end