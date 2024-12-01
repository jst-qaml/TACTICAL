function [lb, ub] = obtainBound(sps_weight_info)
% obtainBound function can determine the lower bound and the upper bound of the suspicious weight in the same
% layer, based on the value of suspicious weights and their positions.
%
% Inputs:
%   sps_weight_info: includes suspiciousness score of a suspicious weight,
%   weight value, layer_idx, right_endpoint_idx, left_endpoint_idx.
% Outputs:
%   lb: lower bound
%   ub: upper bound

[row, ~] = size(sps_weight_info);

lb = zeros(1,row);
ub = zeros(1,row);

% obtain the id of all classes
layers = unique(sps_weight_info(:,3));

% classfy the data and calculate the max and min value
min_layer = min(sps_weight_info(:,3));
max_weight = accumarray(sps_weight_info(:,3) - min_layer + 1, sps_weight_info(:,2), [], @max);
min_weight = accumarray(sps_weight_info(:,3) - min_layer + 1, sps_weight_info(:,2), [], @min);

% temporarily use this method
for i = 1: row

    idx = find(layers == sps_weight_info(i, 3));

    max_w = max_weight(idx, 1);
    min_w = min_weight(idx, 1);

    if min_w < 0
        lb(1,i) = 2 * min_w;
    else
        lb(1,i) = min_w/2;
    end
    if max_w > 0
        ub(1,i) = 2 * max_w;
    else
        ub(1,i) = max_w/2;
    end
end

end