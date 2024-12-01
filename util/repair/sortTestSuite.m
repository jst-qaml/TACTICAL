function sortedTS = sortTestSuite(TS)
% sortTestSuite function sorts the test cases according to their robustness values and
% returns the sorted test suite.
%
% Inputs:
%   TS: the original unsorted test suite
% Outputs:
%   sortedTS: the sorted TS (posTS then negTS)

pos_in = {};
pos_rob = [];
neg_in = {};
neg_rob = [];

% classify the test suite by the robustness value
for i = 1:TS.size
    if TS.tr_rob_set(i,1) > 0
        pos_in{end+1} = TS.tr_in_cell{1,i};
        pos_rob = [pos_rob, TS.tr_rob_set(i,1)];  
    else
        neg_in{end+1} = TS.tr_in_cell{1,i};
        neg_rob = [neg_rob, TS.tr_rob_set(i,1)];
    end
end

% Tests are sorted by decreasing robustness (rob_pt_1 is the highest)
[sorted_pos_rob, sorted_pos_idx] = sort(pos_rob, 'descend');
% reorder the test cases with positive robustness values using the sorted indexes
sorted_pos_in = pos_in(sorted_pos_idx);

% tests are sorted by decreasing robustness (rob_ft_1 is the highest)
[sorted_neg_rob, sorted_neg_idx] = sort(neg_rob, 'descend');
% reorder the test cases with negative robustness values using the sorted indexes
sorted_neg_in = neg_in(sorted_neg_idx);

sortedTS.tr_in_cell = [sorted_pos_in, sorted_neg_in];
sortedTS.tr_rob_set = [sorted_pos_rob, sorted_neg_rob]';
sortedTS.sn = TS.sn;
sortedTS.size = TS.size;
end