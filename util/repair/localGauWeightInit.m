function X = localGauWeightInit(sps_weight_info, pop_size, pert_level)
% localGauWeightInit function can initialize the weights of each candidate
% individual. The initial value of each neural weight (i.e. a suspicious neural weight)
% is sampled from a Gaussian distribution, whose mean and standard deviation are
% computed using all the weights in the same layer with this suspicious neural weight.
%
% Inputs:
%   sps_weight_info: includes suspiciousness score of a suspicious weight,
%   weight value, layer_idx, right_endpoint_idx, left_endpoint_idx.
%   pop_size: the size of population
%   pert_level: perturbation level in [0,1]
% Outputs:
%   X: an initialized population

[row, ~] = size(sps_weight_info);

% initialize X, each row represents a candidate individual
X = repmat(sps_weight_info(:,2)', pop_size, 1);
% add perturbation
for i = 1:pop_size
    for j = 1:row
        sigma = X(i, j) * pert_level;
        % add perturbation to current weight
        X(i, j) = X(i, j) + sigma * randn();
    end
end

% replace the first row as original weight, to guarantee that the worst
% solution will be not lower than that of the original weight
X(1, :) = sps_weight_info(:,2)';
end
