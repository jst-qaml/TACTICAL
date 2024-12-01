function [X, gaussianInfo] = globalGauWeightInit(weight, sps_weight_info, pop_size)
% globalGauWeightInit function can initialize the weights of each candidate
% individual. The initial value of each neural weight (i.e. a suspicious neural weight)
% is sampled from a Gaussian distribution, whose mean and standard deviation are
% computed using all the weights in the same layer with this suspicious neural weight.
%
% Inputs:
%   weight: neural weight of hidden layers and output layer stored as a cell array according to the layer index
%   sps_weight_info: includes suspiciousness score of a suspicious weight,
%   weight value, layer_idx, right_endpoint_idx, left_endpoint_idx.
%   pop_size: the size of population
% Outputs:
%   X: an initialized population
%   gaussianInfo: gaussian distribution information (mean and std) of neural weights in the same layer

[row, ~] = size(sps_weight_info);

% initialize X, each row represents a candidate individual
X = zeros(pop_size, row);

% obtain the layer number of the neural network
layer_num = numel(weight);
gaussianInfo = zeros(layer_num,2);

% calculate the mean and the std of each layer
for li = 1:layer_num
    % the weight info of current layer in the form of mat
    layer_weight = weight{1,li};
    % mean
    gaussianInfo(li,1) = mean(layer_weight,'all');
    % std
    gaussianInfo(li,2) = std(layer_weight(:));
end

% generate X by column
for i = 1:row
    cur_layer_idx = sps_weight_info(i, 3);
    w_column = normrnd(gaussianInfo(cur_layer_idx, 1), gaussianInfo(cur_layer_idx, 2), [pop_size, 1]);
    X(:, i) = w_column;
end

% replace the first row as original weight, to guarantee that the worst
% solution will be not lower than that of the original weight
X(1, :) = sps_weight_info(:,2)';
end
