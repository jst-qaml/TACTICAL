% This script randomly selects 1 mutant from 80% faulty benchmarks to investigate RQ6.
clear;
clc;
close all;
%%
path = '/Users/ldy/git/TACTICAL';
benchList = {
    'ACC2_spec1_simlog', ...
    'AFC2_spec1_simlog', ...
    'WT2_spec1_simlog', ...
    'SC2_spec1_simlog'
    };

fixedMutNum = 1;

for xi = 1:3
    for bmi = 1:4
        benchPath = [path, '/MUT', num2str(xi), '/', benchList{1,bmi}, '/RQ234'];
        dirOutput = dir(fullfile([benchPath, '/'],'*_M_*.mat'));
        filePaths = {dirOutput.folder};
        fileNames = {dirOutput.name};
        mutNum = numel(fileNames);

        RQ6Idx = randperm(mutNum, fixedMutNum);
        % create directories for both parts
        cd(benchPath);
        cd ../../..;
        newFolder = ['RQ6/MUT', num2str(xi), '/', benchList{1,bmi}];
        mkdir(newFolder);
        % copy fixedMutNum files to the RQ5 folder
        for i = 1:numel(RQ6Idx)
            srcFile = fullfile(filePaths{RQ6Idx(i)}, fileNames{RQ6Idx(i)});
            copyfile(srcFile, newFolder);
        end
    end
end

srList = zeros(3,4);

for xi = 1:3
    for bmi = 1:4
        benchPath = [path, '/RQ6/MUT', num2str(xi), '/', benchList{1,bmi}];
        dirOutput = dir(fullfile([benchPath, '/'],'*_M_*.mat'));
        filePaths = {dirOutput.folder};
        fileNames = {dirOutput.name};
        srcFile = fullfile(filePaths{1}, fileNames{1});
        curMut = load(srcFile);
        srList(xi,bmi) = curMut.cur_safety_rate;
    end
end