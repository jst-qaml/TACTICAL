function parSaveSimLog(simlog_file, mut_info, cur_state_set, delta_state, cur_safety_rate, cur_diagInfo_suite)
% parSaveSimLog function can save variables as a file during parallel
% computing.
%
% Inputs:
%   simlog_file: filename
%   the variables need to be saved
% Outputs:

save(simlog_file, 'mut_info', 'cur_state_set', 'delta_state', 'cur_safety_rate', 'cur_diagInfo_suite');
end