function [recallList] = genRecallList(sortedSpsNeuronList, mut_op)
% genRecallList function returns recallList.
%
% Inputs:
%   sortedSpsNeuronList: sorted suspicious neuron list
%   mut_op: mutated neurons
% Outputs:
%   recallList: 0,1/3,1/2,2/3,1

spsNeuNum = size(mut_op,1);
neuNum = size(sortedSpsNeuronList,1);

recallList = zeros(neuNum,1);

for topsi = 1:neuNum
    [detNeuron, ~, ~] = intersect(sortedSpsNeuronList(1:topsi,:), mut_op, 'rows');
    recallList(topsi,1) = size(detNeuron,1)/spsNeuNum;
end

if recallList
end