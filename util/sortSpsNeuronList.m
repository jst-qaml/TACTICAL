function [sortedSpsNeuronList] = sortSpsNeuronList(spsNeuronList, mode)
% sortSpsNeuronList function returns sorted SpsNeuronList.
%
% Inputs:
%   spsNeuronList: suspicious neuron list
%   mode: 'random'/'fixed'
% Outputs:
%   sortedSpsNeuronList: sorted suspicious neuron list

% replace NaN with 0
spsNeuronList(isnan(spsNeuronList)) = 0;

% the index of the column to sort by
colIdx = 3;

% sort the matrix A by colIdx
[sortedSpsNeuronList, sortIdx] = sortrows(spsNeuronList, -colIdx);

if strcmp(mode, 'fixed')
    return;
end

% find unique values in the specified column
uniqueVals = unique(sortedSpsNeuronList(:, colIdx));

% loop through each unique value
for i = 1:length(uniqueVals)
    % find the row indices of the same value
    sameValRows = find(sortedSpsNeuronList(:, colIdx) == uniqueVals(i));

    % if there are multiple rows with the same value, randomize their order
    if numel(sameValRows) > 1
        randOrder = randperm(numel(sameValRows));
        sortedSpsNeuronList(sameValRows, :) = sortedSpsNeuronList(sameValRows(randOrder), :);
    end
end
end



